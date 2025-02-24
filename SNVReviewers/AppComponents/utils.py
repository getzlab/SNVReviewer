import argparse
import pandas as pd
import plotly.graph_objects as go
import numpy as np
import scipy as sp
from statsmodels.stats.multitest import fdrcorrection
import json
import plotly.express as px

# minimum number of rows (genes) to display in the table
n_rows_min = 50
# horizontal buffer for the scatter plots
hor_buffer = 0.01
# vertical gap between gene labels and scatter points
y_gap_annot = 0.01
# buffer for the number of rows in the table
n_rows_buffer = 0.5
# maximum value along the vertical axis for the volcano and Q-Q plots
ymax = 16
# beta confidence interval
ci = 0.95

# Combined dig report properties of significant points
col_sig = 'rgba(255, 0, 0, 1)'
opac_sig = 0.7

# properties of non-significant points
col_nonsig = 'rgba(0, 0, 0, 1)'
opac_nonsig = 0.5
# properties of thinner lines
thk_thin = 0.75
col_thin = 'gray'
typ_thin = 'dash'
# properties of thicker lines
thk_thick = 2.5
col_thick = 'gray'
typ_thick = 'dash'
# properties of error bars
wid_err = 2
thk_err = 0.3
opac_err = 0.3
# properties of the bar plot
col_bar = 'gray'
opac_bar = 0.8

# color of beta confidence area
col_beta = 'rgba(128, 128, 128, 0.2)'

# markers used to indicate the dominant interval set in the QQ plot
markers = {
    'coding': ["Coding", "circle"],
    'promoter': ["Promoter", "square"],
    '5utr': ["5\' UTR", "star"],
    '3utr': ["3\' UTR", "triangle-down"]
}
msize = 7.5

# dropdown options
combined_burden_plot_type = {
    'total': '',
    'sample_wise': 'SAMPLE',
}
combined_result_types = ['coding', 'promoter', '5utr', '3utr']

display_labels_type = {
    'Yes': True,
    'No': False
}

display_bounds_type = {
    'Yes': True,
    'No': False
}

# threshold for log2 fold change (observed vs expected)
logfc_thr = 1
# minimum number of rows (genes) to display in the table
n_rows_min = 50
# horizontal buffer for the scatter plots
hor_buffer = 0.05
# buffer for the number of rows in the table
n_rows_buffer = 0.5
# maximum value along the vertical axis for the volcano and Q-Q plots
ymax = 16

# Coding region properties of significant points
col_sig = 'rgba(255, 0, 0, 1)'
opac_sig = 0.7
# properties of non-significant points
col_nonsig = 'rgba(0, 0, 0, 1)'
opac_nonsig = 0.5
# properties of thinner lines
thk_thin = 0.75
col_thin = 'gray'
typ_thin = 'dash'
# properties of thicker lines
thk_thick = 2.5
col_thick = 'gray'
typ_thick = 'dash'
# properties of error bars
wid_err = 1
thk_err = 0.5
opac_err = 0.3
# properties of the bar plot
col_bar = 'gray'
opac_bar = 0.8
# text for the special case when Sample-wise case does not exist
SPECIAL_TEXT = 'Sample-wise case does not exist for Indels and Indels + Nonsynonymous SNVs!'

# derived parameters
col_err_sig = ','.join(col_sig.split(',')[:-1]) + ', {})'.format(opac_err)
col_err_nonsig = ','.join(col_nonsig.split(',')[:-1]) + ', {})'.format(opac_err)

# dictionaries for coding region the two dropdowns
coding_region_burden_type = {
    'total': 'BURDEN',
    'sample_wise': 'BURDEN_SAMPLE',
}

coding_region_mutation_type = {
    'indels_nonsynonymous_snvs': 'MUT',
    'indels': 'INDEL',
    'nonsynonymous_snvs': 'NONSYN',
    'missense_snvs': 'MIS',
    'nonsense_snvs': 'NONS',
    'truncating_snvs': 'TRUNC',
    'splice_site_snvs': 'SPL',
    'synonymous_snvs': 'SYN',
}

# used for generating the Combined QQ plot
combined_mutation_type = {
    'indels_snvs': 'MUT',
    'indels': 'INDEL',
    'snvs': 'SNV'
}

# same for the combined and coding region
scatterpoint_type = {
    "uniform_p_mid": "unif",
    "p_mid": "recalc"
}


# dictionaries for coding region the two dropdowns
non_coding_region_burden_type = {
    'total': 'BURDEN',
    'sample_wise': 'SAMPLE',
}

def nb_pvalue_greater_midp(k, alpha, p):
    """ Calculate an UPPER TAIL p-value for a negative binomial distribution with a midp correction
    """
    return 0.5 * sp.stats.nbinom.pmf(k, alpha, p) + sp.special.betainc(k + 1, alpha, 1 - p)

def nb_pvalue_lower(k, alpha, p):
    """ Calculate the upper bound for the p-value of a negative binomial distribution.
    """
    return sp.special.betainc(k + 1, alpha, 1 - p)

def nb_pvalue_upper(k, alpha, p):
    """ Calculate the upper bound for the p-value of a negative binomial distribution.
    """
    ind_0 = k == 0
    pvals = np.zeros_like(alpha)
    pvals[ind_0] = sp.stats.nbinom.pmf(k[ind_0], alpha[ind_0], p[ind_0])
    pvals[~ind_0] = sp.special.betainc(k[~ind_0], alpha[~ind_0], 1 - p[~ind_0])
    return pvals

def nb_pvalue_uniform_midp(k, alpha, p):
    """ Calculate the upper tail p-value for negative binomial distribution using uniform approximation and a random draw.
    """
    return np.random.uniform(size=k.shape) * sp.stats.nbinom.pmf(k, alpha, p) + sp.special.betainc(k + 1, alpha, 1 - p)

def coding_region_reformat_numbers(x, format='{:.3E}'):
    """
    Reformat numbers in an array to a specific format
    """
    return [format.format(n) for n in x]


# DIG REPORT Combined Functions
def generate_dig_report_dataframe_combined(
        path_to_dig_results, 
    ):
    """
    Generates the dataframe that contains the data from the DIG report
    """
    
    # Driver gene lists
    cgc_list_path = "gs://getzlab-workflows-reference_files-oa/hg19/dig/cancer_gene_census_2024_06_20.tsv"
    pancan_list_path = "gs://getzlab-workflows-reference_files-oa/hg19/dig/pancanatlas_genes.tsv"

    cgc_list = pd.read_csv(cgc_list_path, sep='\t').to_numpy().flatten()
    pancan_list = pd.read_csv(pancan_list_path, sep='\t').to_numpy().flatten()

    # Output from DIGDriver
    print("this is the path to the dig results: ")
    print(path_to_dig_results)
    df = pd.read_csv(path_to_dig_results, sep='\t')
    # df = df.iloc[:20]
    # Adding indicator of genes being part of the CGC or PanCan list
    df['CGC'] = df['GENE'].isin(cgc_list)
    df['PANCAN'] = df['GENE'].isin(pancan_list)
    muts_ts = list(combined_mutation_type.values())

    if 'EXP_INDEL' in df.columns:

        # Adding new columns for Non-synonymous SNVs + Indels
        df['OBS_MUT'] = df['OBS_NONSYN'] + df['OBS_INDEL']
        df['EXP_MUT'] = df['EXP_NONSYN'] + df['EXP_INDEL']
    else:
        for key in list(combined_mutation_type.keys()):
            if 'indel' in key.lower():
                del combined_mutation_type[key]

    # Computing lower and upper bounds for the p-values
    muts_ts.remove('INDEL')
    muts_ts.remove('MUT')

    for m in muts_ts:
        # total burden
        print("this is the value for m: ", m)
        
        df['PVAL_' + m + '_BURDEN_recalc'] = nb_pvalue_greater_midp(
            df['OBS_' + m],
            df.ALPHA,
            1 / (df.THETA * df['Pi_' + m] + 1)
        )
        df['PVAL_' + m + '_BURDEN_unif'] = nb_pvalue_uniform_midp(
            df['OBS_' + m],
            df.ALPHA,
            1 / (df.THETA * df['Pi_' + m] + 1)
        )
        df['PVAL_' + m + '_BURDEN_lower'] = nb_pvalue_lower(
            df['OBS_' + m],
            df.ALPHA,
            1 / (df.THETA * df['Pi_' + m] + 1)
        )
        df['PVAL_' + m + '_BURDEN_upper'] = nb_pvalue_upper(
            df['OBS_' + m],
            df.ALPHA,
            1 / (df.THETA * df['Pi_' + m] + 1)
        )
        # sample-wise burden
        df['PVAL_' + m + '_BURDEN_SAMPLE_recalc'] = nb_pvalue_greater_midp(
            df['N_SAMP_' + m],
            df.ALPHA,
            1 / (df.THETA * df['Pi_' + m] + 1)
        )
        df['PVAL_' + m + '_BURDEN_SAMPLE_unif'] = nb_pvalue_uniform_midp(
            df['N_SAMP_' + m],
            df.ALPHA,
            1 / (df.THETA * df['Pi_' + m] + 1)
        )
        df['PVAL_' + m + '_BURDEN_SAMPLE_lower'] = nb_pvalue_lower(
            df['N_SAMP_' + m],
            df.ALPHA,
            1 / (df.THETA * df['Pi_' + m] + 1)
        )
        df['PVAL_' + m + '_BURDEN_SAMPLE_upper'] = nb_pvalue_upper(
            df['N_SAMP_' + m],
            df.ALPHA,
            1 / (df.THETA * df['Pi_' + m] + 1)
        )
    if 'EXP_INDEL' in df.columns:
        # total indel burden
        df['PVAL_INDEL_BURDEN_recalc'] = nb_pvalue_greater_midp(
            df.OBS_INDEL,
            df.ALPHA_INDEL,
            1 / (df.THETA_INDEL * df.Pi_INDEL + 1)
        )
        df['PVAL_INDEL_BURDEN_unif'] = nb_pvalue_uniform_midp(
            df.OBS_INDEL,
            df.ALPHA_INDEL,
            1 / (df.THETA_INDEL * df.Pi_INDEL + 1)
        )
        df['PVAL_INDEL_BURDEN_lower'] = nb_pvalue_lower(
            df.OBS_INDEL,
            df.ALPHA_INDEL,
            1 / (df.THETA_INDEL * df.Pi_INDEL + 1)
        )
        df['PVAL_INDEL_BURDEN_upper'] = nb_pvalue_upper(
            df.OBS_INDEL,
            df.ALPHA_INDEL,
            1 / (df.THETA_INDEL * df.Pi_INDEL + 1)
        )
        # p-values for nonsynonymous SNVs + indels
        col_mut = 'PVAL_MUT_BURDEN'
        df[col_mut + '_recalc'] = np.nan
        df[col_mut + '_unif'] = np.nan
        df[col_mut + '_lower'] = np.nan
        df[col_mut + '_upper'] = np.nan

        for idx in df.index:
            df.at[idx, col_mut + '_recalc'] = sp.stats.combine_pvalues(
                [df.at[idx, 'PVAL_NONSYN_BURDEN_recalc'], df.at[idx, 'PVAL_INDEL_BURDEN_recalc']],
                method='fisher')[1]
            df.at[idx, col_mut + '_unif'] = sp.stats.combine_pvalues(
                [df.at[idx, 'PVAL_NONSYN_BURDEN_unif'], df.at[idx, 'PVAL_INDEL_BURDEN_unif']],
                method='fisher')[1]
            df.at[idx, col_mut + '_lower'] = sp.stats.combine_pvalues(
                [df.at[idx, 'PVAL_NONSYN_BURDEN_lower'], df.at[idx, 'PVAL_INDEL_BURDEN_lower']],
                method='fisher')[1]
            df.at[idx, col_mut + '_upper'] = sp.stats.combine_pvalues(
                [df.at[idx, 'PVAL_NONSYN_BURDEN_upper'], df.at[idx, 'PVAL_INDEL_BURDEN_upper']],
                method='fisher')[1]
            
    return df

def plot_qq(x, y, erry, labels, col, opac, mark, name, fig, hi='name+text'):
    fig.add_trace(
        go.Scatter(
            x=x,
            y=y,
            error_y=erry,
            mode='markers',
            marker=dict(color=col, opacity=opac, symbol=mark, size=msize),
            text=labels,
            name=name,
            hoverinfo=hi,
            showlegend=False
        )
    )
    return fig

def get_labels_w_err(x, y, y_lower, y_upper, labels, form="({:.2f}, {:.2f} +{:.2f} / -{:.2f})<br>{}"):
    return [form.format(x[i], y[i], y_upper[i], y_lower[i], l) for i, l in enumerate(labels)]


def get_labels(x, y, labels, form="({:.2f}, {:.2f})<br>{}"):
    return [form.format(x[i], y[i], l) for i, l in enumerate(labels)]

def nb_pvalue_greater_midp(k, alpha, p):
    """ Calculate an UPPER TAIL p-value for a negative binomial distribution
        with a midp correction
    """
    return 0.5 * sp.stats.nbinom.pmf(k, alpha, p) + sp.special.betainc(k+1, alpha, 1-p)

def nb_pvalue_lower(k, alpha, p):
    """ Calculate the upper bound for the p-value of a negative binomial distribution.
    """
    return sp.special.betainc(k+1, alpha, 1-p)

def nb_pvalue_upper(k, alpha, p):
    """ Calculate the upper bound for the p-value of a negative binomial distribution.
    """
    ind_0 = k == 0
    pvals = np.zeros_like(alpha)
    pvals[ind_0] = sp.stats.nbinom.pmf(k[ind_0], alpha[ind_0], p[ind_0])
    pvals[~ind_0] = sp.special.betainc(k[~ind_0], alpha[~ind_0], 1-p[~ind_0])
    return pvals

def nb_pvalue_uniform_midp(k, alpha, p):
    """ Calculate the upper tail p-value for negative binomial distribution using uniform approximation and a random draw.
    """
    return np.random.uniform(size=k.shape) * sp.stats.nbinom.pmf(k, alpha, p) + sp.special.betainc(k+1, alpha, 1-p)

def combined_reformat_numbers(df, cols, form='{:.2E}'):
    """
    Reformat numbers in an array to a specific format
    """
    return df[cols].applymap(lambda x: form.format(x) if not pd.isna(x) else 'NA')

# Make all the markers the same, i.e. don't have the dominant region, make sure to redownload the code from the google drive
def generate_combined_plot_data(df, mut, bur, display_bounds, scatterpoint):
    """
    Given a mutation type and a burden type, generate the data for the volcano plot, Q-Q plot, and table plot
    :param mut: str, mutation type
    :param bur: str, burden type
    :param display_bounds: bool, whether to display the bounds of the p-values
    :param scatterpoint: str, type of p-values to use
    """
    df_comb = df.copy()
    alp = 0.1 # significance level

    col_chosen = mut + ('_SAMPLE' if bur=='SAMPLE' else '') + '_' + scatterpoint
    col_pvals = ['PVAL_' + rt + '_' + col_chosen for rt in combined_result_types]
    col_sizes = ['SIZE_' + rt  for rt in combined_result_types]
    cols_kept = ['GENE', 'CHROM'] + col_pvals + col_sizes + ['PVAL_' + col_chosen] + ['FDR_' + col_chosen] + ['CGC', 'PANCAN']
    
    if display_bounds:
        cols_bound = ['PVAL_' + '_'.join(col_chosen.split('_')[:-1]) + '_' + typ for typ in ['lower', 'upper']]
        df_kept = df_comb[cols_kept + cols_bound].sort_values(by='PVAL_' + col_chosen, ignore_index=True)
        pval_bounds = df_kept[cols_bound].to_numpy()

    else:
        df_kept = df_comb[cols_kept].sort_values(by='PVAL_' + col_chosen, ignore_index=True)
        pval_bounds = None
    df_kept['RANK'] = df_kept.index + 1
    pvals = df_kept['PVAL_' + col_chosen].to_numpy()
    labels = df_kept['GENE'].to_numpy()

    # Determine dominant result type
    test_dom = df_kept[col_pvals].idxmin(axis=1).str.split('_', expand=True)[1].to_numpy()

    # Determine significant points
    ind_sig = df_kept['FDR_' + col_chosen] < alp

    # Dataframe for the plot
    cols_plot = ['RANK', 'GENE', 'CHROM', 'FDR_' + col_chosen, 'PVAL_' + col_chosen] + col_pvals +  col_sizes + ['CGC', 'PANCAN']
    n_rows = max(n_rows_min, int(np.sum(ind_sig) * (1 + n_rows_buffer)))
    df_plot = df_kept.iloc[:n_rows][cols_plot].copy()
    df_plot.rename(columns={
        'PVAL_' + col_chosen: 'PVAL',
        'FDR_' + col_chosen : 'FDR'
    }, inplace=True)
    df_plot.rename(columns={'PVAL_' + rt + '_' + col_chosen: 'PVAL_' + rt for rt in combined_result_types}, inplace=True)

    df_plot[col_sizes] = combined_reformat_numbers(df_plot, col_sizes, form='{:.0f}')
    cols_floats = ['PVAL', 'FDR'] + ['PVAL_' + rt for rt in combined_result_types]
    df_plot[cols_floats] = combined_reformat_numbers(df_plot, cols_floats)

    # Ensure no NaN values in the table
    # df_plot.fillna('NA', inplace=True)
    df_plot[['RANK', 'GENE', 'CGC', 'PANCAN']] = df_plot[['RANK', 'GENE', 'CGC', 'PANCAN']].astype(str)
    # Generate table figure
    headerColor = 'grey'
    rowEvenColor = 'lightgrey'
    rowOddColor = 'white'

    df_plot = df_plot.astype(str)
    for i in range(df_plot.shape[0]):
        if ind_sig[i]:
            df_plot.loc[i, :] = '<b>' + df_plot.loc[i, :] + '</b>'

    # Adding hyperlinks to a Google search for the gene names
    gene_entries = []
    for g in df_plot['GENE']:
        if '<b>' in g:
            g_trimmed = g.split('>')[1].split('<')[0]
        else:
            g_trimmed = g
        gene_entries.append(
            f'<a href="https://www.google.com/search?q={g_trimmed}+gene+cancer" target="_blank">{g}</a>')
    # df_plot['GENE'] = gene_entries
    # df_plot = df_plot[['RANK', 'GENE']]

    table_fig = go.Figure(data=[go.Table(
        header=dict(values=['<b>' + col.replace('_', '<br>') + '</b>' for col in df_plot.columns],
                    line_color='darkslategray',
                    fill_color=headerColor,
                    align=['left'] + ['center'] * (len(df_plot.columns)-1),
                    font=dict(color='white', size=12)
                    ),
        cells=dict(values=[df_plot[col].tolist() for col in df_plot.columns],
                    line_color='darkslategray',
                    fill_color=[[rowOddColor if i % 2 == 0 else rowEvenColor for i in range(df_plot.shape[0])]],
                    align=['left'] + ['center'] * (len(df_plot.columns)-1),
                    font=dict(color='darkslategray', size=11),
                    # format = ['html'] * len(df_plot.columns)  # Enable HTML formatting
                    )
    )])

    return df_kept, pvals, pval_bounds, labels, test_dom, ind_sig, table_fig


# for combined dig report plots
def generate_combined_dig_report_plots(df, mut_key, bur_key, display_bounds_key, display_labels_key, scatterpoint_key):

    # prepare plot data for all combinations of mut_type and burden_type dropdown options
    plot_data = {}
    
    mut_val = combined_mutation_type[mut_key]
    bur_val = combined_burden_plot_type[bur_key]
    display_bounds_val = display_bounds_type[display_bounds_key]
    scatterpoint_val = scatterpoint_type[scatterpoint_key]
    display_labels_val = display_labels_type[display_labels_key]
    text_special = ""

    if not (mut_key in ['indels', 'indels_snvs'] and bur_key == 'sample_wise'):
        _, pvals, pval_bounds, labels, test_dom, ind_sig, table_fig = generate_combined_plot_data(df, mut_val, 
                                                                                                  bur_val, display_bounds_val, 
                                                                                                  scatterpoint_val)

        # Q-Q Plot
        # Scatter plots
        x = -np.log10(np.arange(1, len(pvals) + 1) / (len(pvals) + 1))
        y = -np.log10(pvals)
        ind_capped = y > ymax
        ind_ncapped = np.logical_and(ind_sig, ~ind_capped)
        ind_nonsig = ~ind_sig
        ylim_upper = min(np.max(y), ymax) * (1 + hor_buffer)

        qq_fig = go.Figure()

        for i in range(len(combined_result_types)):

            ind = test_dom == combined_result_types[i]
            marker_i = markers[combined_result_types[i]][1]
            i_nonsig = np.logical_and(ind, ind_nonsig)
            i_ncapped = np.logical_and(ind, ind_ncapped)
            i_capped = np.logical_and(ind, ind_capped)

            labels_capped = labels[i_capped].tolist()
            x_capped = x[i_capped].tolist()
            y_capped = y[i_capped].tolist()
            x_ncapped = x[i_ncapped].tolist()
            y_ncapped = y[i_ncapped].tolist()
            x_nonsig = x[i_nonsig].tolist()
            y_nonsig = y[i_nonsig].tolist()

            # displays the lower and upper bound
            if display_bounds_val:
                y_upper = -np.log10(pval_bounds[:, 0]) - y
                y_lower = y + np.log10(pval_bounds[:, 1])
                y_upper_ncapped = y_upper[i_ncapped].tolist()
                y_lower_ncapped = y_lower[i_ncapped].tolist()
                y_upper_nonsig = y_upper[i_nonsig].tolist()
                y_lower_nonsig = y_lower[i_nonsig].tolist()
                dict_erry_sig = dict(
                    type='data',
                    symmetric=False,
                    array=y_upper_ncapped,
                    arrayminus=y_lower_ncapped,
                    thickness=thk_err,
                    width=wid_err,
                    color=col_err_sig
                )
                dict_erry_nonsig = dict(
                    type='data',
                    symmetric=False,
                    array=y_upper_nonsig,
                    arrayminus=y_lower_nonsig,
                    thickness=thk_err,
                    width=wid_err,
                    color=col_err_nonsig
                )
                labels_capped = get_labels_w_err(x_capped, y_capped, y_lower[i_capped].tolist(),
                                                    y_upper[i_capped].tolist(), labels_capped)
                labels_ncapped = get_labels_w_err(x_ncapped, y_ncapped, y_lower_ncapped,
                                                    y_upper_ncapped, labels[i_ncapped].tolist())
                labels_nonsig = get_labels_w_err(x_nonsig, y_nonsig, y_lower_nonsig, y_upper_nonsig,
                                                labels[i_nonsig].tolist())
                ylim_upper = min(np.max(y_upper + y), ymax) * (1 + hor_buffer)
            else:
                dict_erry_sig, dict_erry_nonsig = None, None
                labels_capped = get_labels(x_capped, y_capped, labels_capped)
                labels_ncapped = get_labels(x_ncapped, y_ncapped, labels[i_ncapped].tolist())
                labels_nonsig = get_labels(x_nonsig, y_nonsig, labels[i_nonsig].tolist())
                ylim_upper = min(np.max(y), ymax) * (1 + hor_buffer)

            qq_fig = plot_qq(x_nonsig, y_nonsig, dict_erry_nonsig, labels_nonsig, col_nonsig, opac_nonsig,
                                marker_i, 'Non-significant', qq_fig)
            qq_fig = plot_qq(x_ncapped, y_ncapped, dict_erry_sig, labels_ncapped, col_sig, opac_sig,
                                marker_i, 'Significant', qq_fig)
            qq_fig = plot_qq(x_capped, [ymax] * sum(i_capped), None, labels_capped, col_sig,
                                opac_sig, marker_i, 'Significant', qq_fig)

        # Dummy points for legend
        for test_i in markers.keys():
            qq_fig.add_trace(
                go.Scatter(
                    y=[None],
                    mode='markers',
                    marker=dict(
                        color='white',
                        symbol=markers[test_i][1],
                        size=msize * 1.25,
                        line=dict(color='black', width=2)
                    ),
                    name=markers[test_i][0]
                )
            )

        if display_labels_val:
            x_capped = x[ind_capped].tolist()
            x_ncapped = x[ind_ncapped].tolist()
            y_capped = [ymax] * sum(ind_capped)
            y_ncapped = y[ind_ncapped].tolist()
            labels_capped = labels[ind_capped].tolist()
            labels_ncapped = labels[ind_ncapped].tolist()
            for (xi, yi, label) in zip(x_capped + x_ncapped, y_capped + y_ncapped,
                                        labels_capped + labels_ncapped):
                qq_fig.add_annotation(
                    x=xi,
                    y=yi - ylim_upper * y_gap_annot,
                    text=label.split('<br>')[-1],
                    showarrow=False,
                    font=dict(color=col_sig),
                    textangle=-90,
                    xanchor="center",
                    yanchor="top"
                )

        # Line plots
        qq_fig.add_trace(
            go.Scatter(
                x=[0, np.max(x)],
                y=[0, np.max(x)],
                mode='lines',
                line=dict(dash=typ_thick, color=col_thick, width=thk_thick),
                showlegend=False,
                hoverinfo='skip'
            )
        )

        if ymax < ylim_upper:
            qq_fig.add_trace(
                go.Scatter(
                    x=[0, np.max(x) * (1 + hor_buffer)],
                    y=[ymax] * 2,
                    mode='lines',
                    line=dict(dash=typ_thin, color=col_thin, width=thk_thin),
                    showlegend=False,
                    hoverinfo='skip'
                )
            )

        # Area of confidence intervals for the identity line
        xi = np.arange(1, len(pvals) + 1)
        clower = -np.log10(sp.stats.beta.ppf((1 - ci) / 2, xi, xi[::-1]))
        cupper = -np.log10(sp.stats.beta.ppf((1 + ci) / 2, xi, xi[::-1]))
        qq_fig.add_trace(go.Scatter(
            x=np.concatenate([x, x[::-1]]).tolist(),  # Combine x values for fill
            y=np.concatenate([clower, cupper[::-1]]).tolist(),  # Combine y values for fill
            fill='toself',
            fillcolor=col_beta,
            line=dict(color='rgba(255,255,255,0)'),  # No line for the filled area
            showlegend=False,
            hoverinfo='skip'
        ))

        # Formatting the figure
        qq_fig.update_layout(
            title='QQ-Plot of P-values:',
            xaxis_title='Expected -Log10(P-value)',
            yaxis_title='Observed -Log10(P-value)',
            xaxis=dict(range=[0, np.max(x) * (1 + hor_buffer)]),
            yaxis=dict(range=[0, ylim_upper]),
            template='plotly_white',
            legend=dict(
                title=dict(
                    text="Dominant Region:",
                ),
                indentation=10
            )
        )

        # Save figures as separate data
        plot_data[f"{mut_key}_{bur_key}_{display_bounds_key}_{scatterpoint_key}_{display_labels_key}"] = {
            'qq': qq_fig.to_dict(),
            'table': table_fig.to_dict(),
            'text': bur_key + ' Mutation Burden of ' + mut_key,
            'textcolor': 'black-text'
        }
    else:
        qq_fig = go.Figure()
        table_fig = go.Figure()
        text_special = SPECIAL_TEXT

        plot_data[f"{mut_key}_{bur_key}_{display_bounds_key}_{scatterpoint_key}_{display_labels_key}"] = {
            'qq': None,
            'table': None,
            'text': SPECIAL_TEXT,
            'textcolor': 'red-text'
        }
    
    return qq_fig, table_fig, text_special

# for generating the coding region plots
def nb_pvalue_greater_midp(k, alpha, p):
    """ 
    Calculate an UPPER TAIL p-value for a negative binomial distribution with a midp correction
    """
    return 0.5 * sp.stats.nbinom.pmf(k, alpha, p) + sp.special.betainc(k + 1, alpha, 1 - p)

def nb_pvalue_lower(k, alpha, p):
    """ 
    Calculate the upper bound for the p-value of a negative binomial distribution.
    """
    return sp.special.betainc(k + 1, alpha, 1 - p)

def nb_pvalue_upper(k, alpha, p):
    """ 
    Calculate the upper bound for the p-value of a negative binomial distribution.
    """
    ind_0 = k == 0
    pvals = np.zeros_like(alpha)
    pvals[ind_0] = sp.stats.nbinom.pmf(k[ind_0], alpha[ind_0], p[ind_0])
    pvals[~ind_0] = sp.special.betainc(k[~ind_0], alpha[~ind_0], 1 - p[~ind_0])

    return pvals

def nb_pvalue_uniform_midp(k, alpha, p):
    """ Calculate the upper tail p-value for negative binomial distribution using uniform approximation and a random draw.
    """
    return np.random.uniform(size=k.shape) * sp.stats.nbinom.pmf(k, alpha, p) + sp.special.betainc(k + 1, alpha, 1 - p)

# DIG REPORT Coding Region functions
def generate_coding_region_plot_data(
        df,
        mut, 
        bur, 
        display_bounds, 
        scatterpoint
    ):
    """
    Given a mutation type and a burden type, generate the data for the volcano plot, Q-Q plot, and table plot
    :param mut: str, mutation type
    :param bur: str, burden type
    :param display_bounds: bool, whether to display the bounds of the p-values
    :param scatterpoint: str, type of p-values to use
    """
    df = df.copy()
    alp = 0.1
    col_pval = 'PVAL_' + mut + '_' + bur
    col_obs = 'OBS' if (bur == 'BURDEN') else 'N_SAMP'
    cols_lfc = [e + '_' + mut for e in [col_obs, 'EXP']]
    ind_keep = df[cols_lfc[1]] > 0

    # subsetting to only those genes for which the expected nuber of mutations is greater than 0
    df_kept = df.loc[ind_keep].copy()
    df_kept['LOGFC_' + mut + '_' + bur] = np.log2(df_kept[cols_lfc[0]] / df_kept[cols_lfc[1]] + 1)
    # df_kept['FDR_' + mut + '_' + bur] = fdrcorrection(df_kept[col_pval])[1]
    df_kept['FDR_' + mut + '_' + bur + '_' + scatterpoint] = fdrcorrection(df_kept[col_pval + '_' + scatterpoint])[
        1]
    if display_bounds:
        df_kept['FDR_' + mut + '_' + bur + '_lower'] = fdrcorrection(df_kept[col_pval + '_lower'])[1]
        df_kept['FDR_' + mut + '_' + bur + '_upper'] = fdrcorrection(df_kept[col_pval + '_upper'])[1]

    df_kept['dNdS_OBS'] = df_kept[col_obs + '_NONSYN'] / df_kept['OBS_SYN']
    df_kept['dNdS_EXP'] = df_kept['EXP_NONSYN'] / df_kept['EXP_SYN']
    df_kept = df_kept.sort_values(by='PVAL_' + mut + '_' + bur + '_' + scatterpoint, ignore_index=True)
    df_kept['RANK'] = df_kept.index + 1

    labels = df_kept.GENE.to_numpy()
    pvals = df_kept['PVAL_' + mut + '_' + bur + '_' + scatterpoint].to_numpy()
    qvals = df_kept['FDR_' + mut + '_' + bur + '_' + scatterpoint].to_numpy()
    logfc = df_kept['LOGFC_' + mut + '_' + bur].to_numpy()
    if display_bounds:
        pval_bounds = df_kept[
            ['PVAL_' + mut + '_' + bur + '_lower', 'PVAL_' + mut + '_' + bur + '_upper']].to_numpy()
    else:
        pval_bounds = None
    logq = -np.log10(qvals)
    if display_bounds:
        logq_bounds = -np.log10(
            df_kept[['FDR_' + mut + '_' + bur + '_lower', 'FDR_' + mut + '_' + bur + '_upper']].to_numpy())
    else:
        logq_bounds = None

    # Determine significant points
    ind_sig = qvals < alp
    ind_lfc = np.abs(logfc) > logfc_thr
    ind_kept = np.logical_and(ind_sig, ind_lfc)

    # Dataframe for the plot
    cols_kept = [
        'RANK', 'GENE', 'CHROM', 'GENE_LENGTH',
        'PVAL_' + mut + '_' + bur + '_' + scatterpoint,
        'FDR_' + mut + '_' + bur + '_' + scatterpoint,
        col_obs + '_' + mut,
        'EXP_' + mut,
        'MU', 'SIGMA', 'dNdS_OBS', 'dNdS_EXP', 'FLAG', 'CGC', 'PANCAN']
    n_rows = max(n_rows_min, int(np.sum(ind_sig) * (1 + n_rows_buffer)))
    df_plot = df_kept.iloc[:n_rows][cols_kept].copy()
    df_plot.rename(columns={
        'GENE_LENGTH': 'LENGTH',
        'PVAL_' + mut + '_' + bur + '_' + scatterpoint: 'PVAL',
        'FDR_' + mut + '_' + bur + '_' + scatterpoint: 'FDR',
        col_obs + '_' + mut: 'OBS',
        'EXP_' + mut: 'EXP'
    }, inplace=True)

    for col in ['PVAL', 'FDR']:
        df_plot[col] = coding_region_reformat_numbers(df_plot[col].to_numpy())
    for col in ['MU', 'SIGMA']:
        df_plot[col] = coding_region_reformat_numbers(df_plot[col].to_numpy(), format='{:.2f}')
    for col in ['dNdS_EXP', 'EXP']:
        df_plot[col] = coding_region_reformat_numbers(df_plot[col].to_numpy(), format='{:.3f}')

    df_plot['OBS'] = df_plot['OBS'].astype(int)
    is_inf = np.logical_or(np.isinf(df_plot.dNdS_OBS.to_numpy()), np.isnan(df_plot.dNdS_OBS.to_numpy()))
    dnds_obs = coding_region_reformat_numbers(df_plot.loc[~is_inf, 'dNdS_OBS'].to_numpy().copy(), format='{:.3f}')
    df_plot['dNdS_OBS'] = df_plot['dNdS_OBS'].astype(str)
    df_plot.loc[~is_inf, 'dNdS_OBS'] = dnds_obs
    df_plot.loc[is_inf, 'dNdS_OBS'] = 'NA'
    is_flagged = df_plot['FLAG'].astype(str).str.title() == 'True'
    df_plot['FLAG'] = is_flagged
    df_plot.loc[is_flagged, 'GENE'] = df_plot['GENE'][is_flagged] + '*'

    # Generate table figure
    headerColor = 'grey'
    rowEvenColor = 'lightgrey'
    rowOddColor = 'white'
    cols_specific = ['PVAL', 'FDR', 'OBS', 'EXP']
    # Making the significant rows bold
    df_plot = df_plot.astype(str)
    for i in range(df_plot.shape[0]):
        if ind_kept[i]:
            df_plot.loc[i, :] = '<b>' + df_plot.loc[i, :].astype(str) + '</b>'
    # Adding hyperlinks to a Google search for the gene names
    gene_entries = []
    for g in df_plot['GENE']:
        if '<b>' in g:
            g_trimmed = g.split('>')[1].split('<')[0]
        else:
            g_trimmed = g
        gene_entries.append(
            f'<a href="https://www.google.com/search?q={g_trimmed}+gene+cancer" target="_blank">{g}</a>')
    # df_plot['GENE'] = gene_entries
    # df_plot = df_plot[['RANK', 'GENE', 'OBS']]

    table_fig = go.Figure(data=[go.Table(
        header=dict(values=['<b>' + col + '</b>' for col in df_plot.columns],
                    line_color='darkslategray',
                    fill_color=headerColor,
                    align=['left'] + ['center'] * (len(df_plot.columns) - 1),
                    font=dict(color='white', size=12)
                    ),
        cells=dict(values=[df_plot[col].tolist() for col in df_plot.columns],
                    line_color='darkslategray',
                    fill_color=[[rowOddColor if i % 2 == 0 else rowEvenColor for i in range(df_plot.shape[0])]],
                    align=['left'] + ['center'] * (len(df_plot.columns) - 1),
                    font=dict(color='darkslategray', size=11),
                    )
    )
    ])
    table_fig.update_layout(
        annotations=[
            dict(
                text="*FLAG=True: At least one kilobase-scale region overlapped by gene is <50% uniquely mappable or in the top 99.99th percentile of mutation rate.",
                x=0,
                y=-0.15,
                xref="paper",
                yref="paper",
                showarrow=False,
                align="left",
                valign="top",
                font=dict(size=12)
            )
        ]
    )

    return df_kept, pvals, pval_bounds, logfc, logq, logq_bounds, labels, ind_kept, table_fig

def generate_coding_region_dig_dataframe(
        path_to_dig_results, 
        # dir_output, 
        # cgc_list_path, 
        # pancan_list_path, 
        # prefix_output=None, 
        alp=0.1
    ):
    """ 
    
    """
    # Driver gene lists
    cgc_list_path = "gs://getzlab-workflows-reference_files-oa/hg19/dig/cancer_gene_census_2024_06_20.tsv"
    pancan_list_path = "gs://getzlab-workflows-reference_files-oa/hg19/dig/pancanatlas_genes.tsv"

    # Driver gene lists
    cgc_list = pd.read_csv(cgc_list_path, sep='\t').to_numpy().flatten()
    pancan_list = pd.read_csv(pancan_list_path, sep='\t').to_numpy().flatten()
    # Output from DIGDriver
    df = pd.read_csv(path_to_dig_results, sep='\t')
    # df = df.iloc[:20]
    # Adding indicator of genes being part of the CGC or PanCan list
    df['CGC'] = df['GENE'].isin(cgc_list)
    df['PANCAN'] = df['GENE'].isin(pancan_list)
    muts_ts = list(coding_region_mutation_type.values())

    if 'EXP_INDEL' in df.columns:
        # Adding new columns for Non-synonymous SNVs + Indels
        df['OBS_MUT'] = df['OBS_NONSYN'] + df['OBS_INDEL']
        df['EXP_MUT'] = df['EXP_NONSYN'] + df['EXP_INDEL']
    else:
        for key in list(coding_region_mutation_type.keys()):
            if 'indel' in key.lower():
                del coding_region_mutation_type[key]
    # Computing lower and upper bounds for the p-values
    muts_ts.remove('INDEL')
    muts_ts.remove('MUT')

    for m in muts_ts:
        # total burden
        df['PVAL_' + m + '_BURDEN_recalc'] = nb_pvalue_greater_midp(
            df['OBS_' + m],
            df.ALPHA,
            1 / (df.THETA * df['Pi_' + m] + 1)
        )
        df['PVAL_' + m + '_BURDEN_unif'] = nb_pvalue_uniform_midp(
            df['OBS_' + m],
            df.ALPHA,
            1 / (df.THETA * df['Pi_' + m] + 1)
        )
        df['PVAL_' + m + '_BURDEN_lower'] = nb_pvalue_lower(
            df['OBS_' + m],
            df.ALPHA,
            1 / (df.THETA * df['Pi_' + m] + 1)
        )
        df['PVAL_' + m + '_BURDEN_upper'] = nb_pvalue_upper(
            df['OBS_' + m],
            df.ALPHA,
            1 / (df.THETA * df['Pi_' + m] + 1)
        )
        # sample-wise burden
        df['PVAL_' + m + '_BURDEN_SAMPLE_recalc'] = nb_pvalue_greater_midp(
            df['N_SAMP_' + m],
            df.ALPHA,
            1 / (df.THETA * df['Pi_' + m] + 1)
        )
        df['PVAL_' + m + '_BURDEN_SAMPLE_unif'] = nb_pvalue_uniform_midp(
            df['N_SAMP_' + m],
            df.ALPHA,
            1 / (df.THETA * df['Pi_' + m] + 1)
        )
        df['PVAL_' + m + '_BURDEN_SAMPLE_lower'] = nb_pvalue_lower(
            df['N_SAMP_' + m],
            df.ALPHA,
            1 / (df.THETA * df['Pi_' + m] + 1)
        )
        df['PVAL_' + m + '_BURDEN_SAMPLE_upper'] = nb_pvalue_upper(
            df['N_SAMP_' + m],
            df.ALPHA,
            1 / (df.THETA * df['Pi_' + m] + 1)
        )
    if 'EXP_INDEL' in df.columns:
        # total indel burden
        df['PVAL_INDEL_BURDEN_recalc'] = nb_pvalue_greater_midp(
            df.OBS_INDEL,
            df.ALPHA_INDEL,
            1 / (df.THETA_INDEL * df.Pi_INDEL + 1)
        )
        df['PVAL_INDEL_BURDEN_unif'] = nb_pvalue_uniform_midp(
            df.OBS_INDEL,
            df.ALPHA_INDEL,
            1 / (df.THETA_INDEL * df.Pi_INDEL + 1)
        )
        df['PVAL_INDEL_BURDEN_lower'] = nb_pvalue_lower(
            df.OBS_INDEL,
            df.ALPHA_INDEL,
            1 / (df.THETA_INDEL * df.Pi_INDEL + 1)
        )
        df['PVAL_INDEL_BURDEN_upper'] = nb_pvalue_upper(
            df.OBS_INDEL,
            df.ALPHA_INDEL,
            1 / (df.THETA_INDEL * df.Pi_INDEL + 1)
        )
        # p-values for nonsynonymous SNVs + indels
        col_mut = 'PVAL_MUT_BURDEN'
        df[col_mut + '_recalc'] = np.nan
        df[col_mut + '_unif'] = np.nan
        df[col_mut + '_lower'] = np.nan
        df[col_mut + '_upper'] = np.nan
        for idx in df.index:
            df.at[idx, col_mut + '_recalc'] = sp.stats.combine_pvalues(
                [df.at[idx, 'PVAL_NONSYN_BURDEN_recalc'], df.at[idx, 'PVAL_INDEL_BURDEN_recalc']],
                method='fisher')[1]
            df.at[idx, col_mut + '_unif'] = sp.stats.combine_pvalues(
                [df.at[idx, 'PVAL_NONSYN_BURDEN_unif'], df.at[idx, 'PVAL_INDEL_BURDEN_unif']],
                method='fisher')[1]
            df.at[idx, col_mut + '_lower'] = sp.stats.combine_pvalues(
                [df.at[idx, 'PVAL_NONSYN_BURDEN_lower'], df.at[idx, 'PVAL_INDEL_BURDEN_lower']],
                method='fisher')[1]
            df.at[idx, col_mut + '_upper'] = sp.stats.combine_pvalues(
                [df.at[idx, 'PVAL_NONSYN_BURDEN_upper'], df.at[idx, 'PVAL_INDEL_BURDEN_upper']],
                method='fisher')[1]


    return df
    
def generate_coding_region_report(
        df,
        mut_key, 
        bur_key, 
        display_bounds_key, 
        scatterpoint_key,
        alp=0.1
    ):
    # prepare plot data for all combinations of mut_type and burden_type dropdown options
    plot_data = {}
    # alp = 0.1
    df = df.copy()
    text_special = ""
    mut_val = coding_region_mutation_type[mut_key]
    bur_val = coding_region_burden_type[bur_key]
    display_bounds_val = display_bounds_type[display_bounds_key]
    scatterpoint_val = scatterpoint_type[scatterpoint_key]

    if not (mut_key in ['indels', 'indels_nonsynonymous_snvs'] and bur_key == 'sample_wise'):
        df_kept, pvals, pval_bounds, logfc, logq, logq_bounds, labels, ind_kept, table_fig = generate_coding_region_plot_data(df, 
                                                                                                                              mut_val, 
                                                                                                                              bur_val, 
                                                                                                                              display_bounds_val, 
                                                                                                                              scatterpoint_val)

        # Volcano Plot
        # Scatter plots
        ind_capped = logq > ymax
        ind_ncapped = np.logical_and(ind_kept, ~ind_capped)
        labels_capped = labels[ind_capped].tolist()
        logfc_capped = logfc[ind_capped].tolist()
        logq_capped = logq[ind_capped].tolist()

        if display_bounds_val:
            logq_upper = logq_bounds[:, 0] - logq
            logq_lower = logq - logq_bounds[:, 1]
            dict_erry_sig = dict(
                type='data',
                symmetric=False,
                array=np.round(logq_upper[ind_ncapped], 3).tolist(),
                arrayminus=np.round(logq_lower[ind_ncapped], 3).tolist(),
                thickness=thk_err,
                width=wid_err,
                color=col_err_sig
            )
            dict_erry_nonsig = dict(
                type='data',
                symmetric=False,
                array=np.round(logq_upper[~ind_kept], 3).tolist(),
                arrayminus=np.round(logq_lower[~ind_kept], 3).tolist(),
                thickness=thk_err,
                width=wid_err,
                color=col_err_nonsig
            )
            logq_upper_capped = logq_upper[ind_capped].tolist()
            logq_lower_capped = logq_lower[ind_capped].tolist()
            labels_capped = [(
                f"({logfc_capped[i]:.3f}, {logq_capped[i]:.3f} +{logq_upper_capped[i]:.3f} / -{logq_lower_capped[i]:.3f})<br>{l}")
                for i, l in enumerate(labels_capped)]
            ylim_upper = min(np.max(logq_upper + logq), ymax) * (1 + hor_buffer)
        else:
            dict_erry_sig, dict_erry_nonsig = None, None
            labels_capped = [(
                f"({logfc_capped[i]:.3f}, {logq_capped[i]:.3f})<br>{l}")
                for i, l in enumerate(labels_capped)]
            ylim_upper = min(np.max(logq), ymax) * (1 + hor_buffer)

        volcano_fig = go.Figure()

        volcano_fig.add_trace(
            go.Scatter(
                x=logfc[~ind_kept].tolist(),
                y=logq[~ind_kept].tolist(),
                error_y=dict_erry_nonsig,
                mode='markers',
                marker=dict(color=col_nonsig, opacity=opac_nonsig),
                text=labels[~ind_kept].tolist(),
                name='Non-significant',
                showlegend=False,
                xhoverformat='.3f',
                yhoverformat='.3f'
            )
        )

        volcano_fig.add_trace(
            go.Scatter(
                x=logfc[ind_ncapped].tolist(),
                y=logq[ind_ncapped].tolist(),
                error_y=dict_erry_sig,
                mode='markers',
                marker=dict(color=col_sig, opacity=opac_sig),
                text=labels[ind_ncapped].tolist(),
                name='Significant',
                showlegend=False,
                xhoverformat='.3f',
                yhoverformat='.3f'
            )
        )

        volcano_fig.add_trace(
            go.Scatter(
                x=logfc[ind_capped].tolist(),
                y=[ymax] * sum(ind_capped),
                mode='markers',
                marker=dict(color=col_sig, opacity=opac_sig),
                text=labels_capped,
                name='Significant',
                showlegend=False,
                hoverinfo='name+text',
            )
        )

        # Line plots
        volcano_fig.add_trace(
            go.Scatter(
                x=[1, 1],
                y=[0, ylim_upper],
                mode='lines',
                line=dict(dash=typ_thin, color=col_thin, width=thk_thin),
                showlegend=False,
                hoverinfo='skip'
            )
        )
        
        volcano_fig.add_trace(
            go.Scatter(
                x=[0, np.max(logfc) * (1 + hor_buffer)],
                y=[-np.log10(alp), -np.log10(alp)],
                mode='lines',
                line=dict(dash=typ_thick, color=col_thick, width=thk_thick),
                showlegend=False,
                hoverinfo='skip'
            )
        )
        if ylim_upper > ymax:
            volcano_fig.add_trace(
                go.Scatter(
                    x=[0, np.max(logfc) * (1 + hor_buffer)],
                    y=[ymax] * 2,
                    mode='lines',
                    line=dict(dash=typ_thin, color=col_thin, width=thk_thin),
                    showlegend=False,
                    hoverinfo='skip'
                )
            )

        # Figure formatting
        volcano_fig.update_layout(
            title='Observed/Expected counts vs. False Discovery Rate:',
            xaxis_title='Log2(Observed/Expected + 1)',
            yaxis_title='-Log10(FDR)',
            xaxis=dict(range=[0, np.max(logfc) * (1 + hor_buffer)]),
            yaxis=dict(range=[0, ylim_upper]),
            template='plotly_white'
        )

        # Q-Q Plot
        # Scatter plots
        x = -np.log10(np.arange(1, len(pvals) + 1) / (len(pvals) + 1))
        y = -np.log10(pvals)
        ind_capped = y > ymax
        ind_ncapped = np.logical_and(ind_kept, ~ind_capped)
        labels_capped = labels[ind_capped].tolist()
        x_capped = x[ind_capped].tolist()
        y_capped = y[ind_capped].tolist()

        # adds the upper and lower bounds to the 
        if display_bounds_val:
            y_upper = -np.log10(pval_bounds[:, 0]) - y
            y_lower = y + np.log10(pval_bounds[:, 1])
            dict_erry_sig = dict(
                type='data',
                symmetric=False,
                array=np.round(y_upper[ind_ncapped], 3).tolist(),
                arrayminus=np.round(y_lower[ind_ncapped], 3).tolist(),
                thickness=thk_err,
                width=wid_err,
                color=col_err_sig
            )
            dict_erry_nonsig = dict(
                type='data',
                symmetric=False,
                array=np.round(y_upper[~ind_kept], 3).tolist(),
                arrayminus=np.round(y_lower[~ind_kept], 3).tolist(),
                thickness=thk_err,
                width=wid_err,
                color=col_err_nonsig
            )
            y_upper_capped = y_upper[ind_capped].tolist()
            y_lower_capped = y_lower[ind_capped].tolist()
            labels_capped = [(
                f"({x_capped[i]:.3f}, {y_capped[i]:.3f} +{y_upper_capped[i]:.3f} / -{y_lower_capped[i]:.3f})<br>{l}")
                for i, l in enumerate(labels_capped)]
            ylim_upper = min(np.max(y_upper + y), ymax) * (1 + hor_buffer)
        else:
            dict_erry_sig, dict_erry_nonsig = None, None
            labels_capped = [(
                f"({x_capped[i]:.3f}, {y_capped[i]:.3f})<br>{l}")
                for i, l in enumerate(labels_capped)]
            ylim_upper = min(np.max(y), ymax) * (1 + hor_buffer)

        xi = np.arange(1, len(pvals) + 1)
        clower = -np.log10(sp.stats.beta.ppf((1 - ci) / 2, xi, xi[::-1]))
        cupper = -np.log10(sp.stats.beta.ppf((1 + ci) / 2, xi, xi[::-1]))

        qq_fig = go.Figure()
        qq_fig.add_trace(
            go.Scatter(
                x=x[~ind_kept].tolist(),
                y=y[~ind_kept].tolist(),
                error_y=dict_erry_nonsig,
                mode='markers',
                marker=dict(color=col_nonsig, opacity=opac_nonsig),
                text=labels[~ind_kept].tolist(),
                name='Non-significant',
                showlegend=False,
                xhoverformat='.3f',
                yhoverformat='.3f'
            )
        )

        qq_fig.add_trace(
            go.Scatter(
                x=x[ind_ncapped].tolist(),
                y=y[ind_ncapped].tolist(),
                error_y=dict_erry_sig,
                mode='markers',
                marker=dict(color=col_sig, opacity=opac_sig),
                text=labels[ind_ncapped].tolist(),
                name='Significant',
                showlegend=False,
                xhoverformat='.3f',
                yhoverformat='.3f'
            )
        )

        qq_fig.add_trace(
            go.Scatter(
                x=x_capped,
                y=[ymax] * sum(ind_capped),
                mode='markers',
                marker=dict(color=col_sig, opacity=opac_sig),
                text=labels_capped,
                name='Significant',
                showlegend=False,
                hoverinfo='name+text',
            )
        )

        # Line plots
        qq_fig.add_trace(
            go.Scatter(
                x=[0, np.max(x)],
                y=[0, np.max(x)],
                mode='lines',
                line=dict(dash=typ_thick, color=col_thick, width=thk_thick),
                showlegend=False,
                hoverinfo='skip'
            )
        )
        if ylim_upper > ymax:
            qq_fig.add_trace(
                go.Scatter(
                    x=[0, np.max(x) * (1 + hor_buffer)],
                    y=[ymax] * 2,
                    mode='lines',
                    line=dict(dash=typ_thin, color=col_thin, width=thk_thin),
                    showlegend=False,
                    hoverinfo='skip'
                )
            )

        # Area of confidence intervals for the identity line
        qq_fig.add_trace(go.Scatter(
            x=np.concatenate([x, x[::-1]]).tolist(),  # Combine x values for fill
            y=np.concatenate([clower, cupper[::-1]]).tolist(),  # Combine y values for fill
            fill='toself',
            fillcolor=col_beta,
            line=dict(color='rgba(255,255,255,0)'),  # No line for the filled area
            showlegend=False,
            hoverinfo='skip'
        ))

        # Figure formatting
        qq_fig.update_layout(
            title='QQ-Plot of P-values:',
            xaxis_title='Expected -Log10(P-value)',
            yaxis_title='Observed -Log10(P-value)',
            xaxis=dict(range=[0, np.max(x) * (1 + hor_buffer)]),
            yaxis=dict(range=[0, ylim_upper]),
            template='plotly_white'
        )

        # dNdS Plot
        # Scatter plots
        ind_isna = np.logical_or(df_kept['dNdS_EXP'].isna(), df_kept['dNdS_OBS'].isna())
        dnds_obs = df_kept['dNdS_OBS'][~ind_isna].to_numpy()
        dnds_exp = df_kept['dNdS_EXP'][~ind_isna].to_numpy()
        dnds_labels = df_kept['GENE'][~ind_isna].to_numpy()
        xmax = np.max(dnds_exp) * (1 + hor_buffer)
        ind_psel = dnds_obs > dnds_exp
        dnds_fig = go.Figure()
        dnds_fig.add_trace(
            go.Scatter(
                x=dnds_exp[~ind_psel].tolist(),
                y=dnds_obs[~ind_psel].tolist(),
                mode='markers',
                marker=dict(color=col_nonsig, opacity=opac_nonsig),
                text=dnds_labels[~ind_psel].tolist(),
                name='Lower than expected',
                showlegend=False
            )
        )
        dnds_fig.add_trace(
            go.Scatter(
                x=dnds_exp[ind_psel].tolist(),
                y=dnds_obs[ind_psel].tolist(),
                mode='markers',
                marker=dict(color=col_sig, opacity=opac_sig),
                text=dnds_labels[ind_psel].tolist(),
                name='Higher than expected',
                showlegend=False
            )
        )

        # Line plot
        dnds_fig.add_trace(
            go.Scatter(
                x=[0, xmax],
                y=[0, xmax],
                mode='lines',
                line=dict(dash=typ_thick, color=col_thick, width=thk_thick),
                showlegend=False,
                hoverinfo='skip'
            )
        )

        # Figure formatting
        dnds_fig.update_layout(
            title='Observed vs Expected dNdS ratios:',
            xaxis_title='Expected dNdS',
            yaxis_title='Observed dNdS',
            xaxis=dict(range=[0, xmax]),
            yaxis=dict(range=[0, np.max(dnds_obs) * (1 + hor_buffer)]),
            template='plotly_white'
        )

        # Save figures as separate data
        plot_data[f"{mut_key}_{bur_key}_{display_bounds_key}_{scatterpoint_key}"] = {
            'volcano': volcano_fig.to_dict(),
            'qq': qq_fig.to_dict(),
            'dnds': dnds_fig.to_dict(),
            'table': table_fig.to_dict(),
            'text': bur_key + ' Mutation Burden of ' + mut_key,
            'textcolor': 'black-text'
        }
    else:
        text_special = SPECIAL_TEXT
        plot_data[f"{mut_key}_{bur_key}_{display_bounds_key}_{scatterpoint_key}"] = {
            'volcano': None,
            'qq': None,
            'dnds': None,
            'table': None,
            'text': SPECIAL_TEXT    ,
            'textcolor': 'red-text'
        }

    # generate histograms for MU and SIGMA
    fig_mu = go.Figure(data=[go.Histogram(
        x=df_kept['MU'],
        opacity=opac_bar,
        marker_color=col_bar
    )])
    fig_mu.update_layout(
        title='Mean of GP model:',
        xaxis_title='MU (mutations per kilobase)',
        yaxis_title='Number of genes',
        template='plotly_white',
        yaxis_type='log'
    )

    fig_sigma = go.Figure(data=[go.Histogram(
        x=df_kept['SIGMA'],
        opacity=opac_bar,
        marker_color=col_bar
    )])
    fig_sigma.update_layout(
        title='Standard deviation of GP model:',
        xaxis_title='SIGMA (mutations per kilobase)',
        yaxis_title='Number of genes',
        template='plotly_white',
        yaxis_type='log'
    )

    return df_kept, volcano_fig, qq_fig, fig_mu, fig_sigma, dnds_fig, table_fig, text_special

# 3 prime utr report results
def generate_dig_non_coding_region_dataframe(
                                path_to_dig_results, 
                                    ):
    # Driver gene lists
    cgc_list_path = "gs://getzlab-workflows-reference_files-oa/hg19/dig/cancer_gene_census_2024_06_20.tsv"
    pancan_list_path = "gs://getzlab-workflows-reference_files-oa/hg19/dig/pancanatlas_genes.tsv"

    cgc_list = pd.read_csv(cgc_list_path, sep='\t').to_numpy().flatten()
    pancan_list = pd.read_csv(pancan_list_path, sep='\t').to_numpy().flatten()
    # Output from DIGDriver
    df = pd.read_csv(path_to_dig_results, sep='\t')
    # Extract gene name and Ensembl ID
    df['GENE'] = df.ELT.str.split('::', expand=True)[2]
    df['ENSEMBL_ID'] = df.ELT.str.split('::', expand=True)[3]
    # df = df.iloc[:20]
    # Adding indicator of genes being part of the CGC or PanCan list
    df['CGC'] = df['GENE'].isin(cgc_list)
    df['PANCAN'] = df['GENE'].isin(pancan_list)
    if 'EXP_INDEL' in df.columns:
        # Adding new columns for Non-synonymous SNVs + Indels
        df['OBS_MUT'] = df['OBS_SNV'] + df['OBS_INDEL']
        df['EXP_MUT'] = df['EXP_SNV'] + df['EXP_INDEL']
        # Computing lower and upper bounds for the p-values
        pfxs_obs = ['SNV', 'INDEL', 'SAMPLES', 'MUT']
        pfxs_pval = ['SNV', 'INDEL', 'SAMPLE', 'MUT']
        pfxs_pi = ['SUM', 'INDEL', 'SUM', 'MUT']
        pfxs_at = ['', '_INDEL', '', '', '']
    else:
        for key in list(combined_mutation_type.keys()):
            if 'indel' in key.lower():
                del combined_mutation_type[key]
        pfxs_obs = ['SNV', 'SAMPLES']
        pfxs_pval = ['SNV', 'SAMPLE']
        pfxs_pi = ['SUM', 'SUM']
        pfxs_at = ['', '']

    for i in range(len(pfxs_obs)):
        if pfxs_obs[i] == 'MUT':
            col_i = 'PVAL_' + pfxs_pval[i] + '_BURDEN'
            df[col_i + '_recalc'] = np.nan
            df[col_i + '_unif'] = np.nan
            df[col_i + '_lower'] = np.nan
            df[col_i + '_upper'] = np.nan
            for idx in df.index:
                df.at[idx, col_i + '_recalc'] = sp.stats.combine_pvalues(
                    [df.at[idx, 'PVAL_SNV_BURDEN_recalc'], df.at[idx, 'PVAL_INDEL_BURDEN_recalc']],
                    method='fisher')[1]
                df.at[idx, col_i + '_unif'] = sp.stats.combine_pvalues(
                    [df.at[idx, 'PVAL_SNV_BURDEN_unif'], df.at[idx, 'PVAL_INDEL_BURDEN_unif']],
                    method='fisher')[1]
                df.at[idx, col_i + '_lower'] = sp.stats.combine_pvalues(
                    [df.at[idx, 'PVAL_SNV_BURDEN_lower'], df.at[idx, 'PVAL_INDEL_BURDEN_lower']],
                    method='fisher')[1]
                df.at[idx, col_i + '_upper'] = sp.stats.combine_pvalues(
                    [df.at[idx, 'PVAL_SNV_BURDEN_upper'], df.at[idx, 'PVAL_INDEL_BURDEN_upper']],
                    method='fisher')[1]

        else:
            df['PVAL_' + pfxs_pval[i] + '_BURDEN_recalc'] = nb_pvalue_greater_midp(
                df['OBS_' + pfxs_obs[i]],
                df['ALPHA' + pfxs_at[i]],
                1 / (df['THETA' + pfxs_at[i]] * df['Pi_' + pfxs_pi[i]] + 1)
            )
            df['PVAL_' + pfxs_pval[i] + '_BURDEN_unif'] = nb_pvalue_uniform_midp(
                df['OBS_' + pfxs_obs[i]],
                df['ALPHA' + pfxs_at[i]],
                1 / (df['THETA' + pfxs_at[i]] * df['Pi_' + pfxs_pi[i]] + 1)
            )
            df['PVAL_' + pfxs_pval[i] + '_BURDEN_lower'] = nb_pvalue_lower(
                df['OBS_' + pfxs_obs[i]],
                df['ALPHA' + pfxs_at[i]],
                1 / (df['THETA' + pfxs_at[i]] * df['Pi_' + pfxs_pi[i]] + 1)
            )
            df['PVAL_' + pfxs_pval[i] + '_BURDEN_upper'] = nb_pvalue_upper(
                df['OBS_' + pfxs_obs[i]],
                df['ALPHA' + pfxs_at[i]],
                1 / (df['THETA' + pfxs_at[i]] * df['Pi_' + pfxs_pi[i]] + 1)
            )

    return df

def generate_non_coding_region_plot_data(df, mut, bur, display_bounds, scatterpoint, alp=0.1):
    """
    Given a mutation type and a burden type, generate the data for the volcano plot, Q-Q plot, and table plot
    :param mut: str, mutation type
    :param bur: str, burden type
    :param display_bounds: bool, whether to display the bounds of the p-values
    :param scatterpoint: str, type of p-values to use
    """
    col_obs = 'OBS_' + mut if (bur == 'BURDEN') else 'OBS_SAMPLES'
    col_exp = 'EXP_' + mut
    col_pval = 'PVAL_' + (mut + '_' if bur == 'BURDEN' else '') + bur

    # subsetting to only those genes for which the expected nuber of mutations is greater than 0
    ind_keep = df[col_exp] > 0
    df_kept = df.loc[ind_keep].copy()

    df_kept['LOGFC_' + mut + '_' + bur] = np.log2(df_kept[col_obs] / df_kept[col_exp] + 1)
    # df_kept['FDR_' + mut + '_' + bur] = fdrcorrection(df_kept[col_pval])[1]
    df_kept['FDR_' + mut + '_' + bur + '_' + scatterpoint] = fdrcorrection(df_kept[col_pval + '_' + scatterpoint])[1]

    if display_bounds:
        df_kept['FDR_' + mut + '_' + bur + '_lower'] = fdrcorrection(df_kept[col_pval + '_lower'])[1]
        df_kept['FDR_' + mut + '_' + bur + '_upper'] = fdrcorrection(df_kept[col_pval + '_upper'])[1]
    df_kept = df_kept.sort_values(by=col_pval + '_' + scatterpoint, ignore_index=True)
    df_kept['RANK'] = df_kept.index + 1

    labels = df_kept.GENE.to_numpy()
    logfc = df_kept['LOGFC_' + mut + '_' + bur].to_numpy()
    pvals = df_kept[col_pval + '_' + scatterpoint].to_numpy()
    if display_bounds:
        pval_bounds = df_kept[[col_pval + '_lower', col_pval + '_upper']].to_numpy()
    else:
        pval_bounds = None
    qvals = df_kept['FDR_' + mut + '_' + bur + '_' + scatterpoint].to_numpy()
    logq = -np.log10(qvals)
    if display_bounds:
        logq_bounds = -np.log10(
        df_kept[['FDR_' + mut + '_' + bur + '_lower', 'FDR_' + mut + '_' + bur + '_upper']].to_numpy())
    else:
        logq_bounds = None

    # Determine significant points
    ind_sig = qvals < alp
    ind_lfc = np.abs(logfc) > logfc_thr
    ind_kept = np.logical_and(ind_sig, ind_lfc)

    # Dataframe for the plot
    cols_kept = [
        'RANK', 'GENE', 'ELT_SIZE',
        col_pval + '_' + scatterpoint,
        'FDR_' + mut + '_' + bur + '_' + scatterpoint,
        col_obs,
        col_exp,
        'MU', 'SIGMA',
        'FLAG', 'CGC', 'PANCAN']
    n_rows = max(n_rows_min, int(np.sum(ind_sig) * (1 + n_rows_buffer)))
    df_plot = df_kept.iloc[:n_rows][cols_kept].copy()
    df_plot.rename(columns={
        'ELT_SIZE': 'SIZE',
        col_pval + '_' + scatterpoint: 'PVAL',
        'FDR_' + mut + '_' + bur + '_' + scatterpoint: 'FDR',
        col_obs: 'OBS',
        col_exp: 'EXP'
    }, inplace=True)

    for col in ['PVAL', 'FDR']:
        df_plot[col] = coding_region_reformat_numbers(df_plot[col].to_numpy())
    for col in ['MU', 'SIGMA']:
        df_plot[col] = coding_region_reformat_numbers(df_plot[col].to_numpy(), format='{:.2f}')
    for col in ['EXP']:
        df_plot[col] = coding_region_reformat_numbers(df_plot[col].to_numpy(), format='{:.3f}')

    is_flagged = df_plot['FLAG'].astype(str).str.title() == 'True'
    df_plot['FLAG'] = is_flagged
    df_plot.loc[is_flagged, 'GENE'] = df_plot['GENE'][is_flagged] + '*'

    # Generate table figure
    headerColor = 'grey'
    rowEvenColor = 'lightgrey'
    rowOddColor = 'white'
    cols_specific = ['PVAL', 'FDR', 'OBS', 'EXP']
    # Making the significant rows bold
    df_plot = df_plot.astype(str)
    for i in range(df_plot.shape[0]):
        if ind_kept[i]:
            df_plot.loc[i, :] = '<b>' + df_plot.loc[i, :].astype(str) + '</b>'
    # Adding hyperlinks to a Google search for the gene names
    gene_entries = []
    for g in df_plot['GENE']:
        if '<b>' in g:
            g_trimmed = g.split('>')[1].split('<')[0]
        else:
            g_trimmed = g
        gene_entries.append(
            f'<a href="https://www.google.com/search?q={g_trimmed}+gene+cancer" target="_blank">{g}</a>')
    # df_plot['GENE'] = gene_entries
    # df_plot = df_plot[['RANK', 'GENE', 'OBS']]

    table_fig = go.Figure(data=[go.Table(
        header=dict(values=['<b>' + col + '</b>' for col in df_plot.columns],
                    line_color='darkslategray',
                    fill_color=headerColor,
                    align=['left'] + ['center'] * (len(df_plot.columns)-1),
                    font=dict(color='white', size=12)
                    ),
        cells=dict(values=[df_plot[col].tolist() for col in df_plot.columns],
                    line_color='darkslategray',
                    fill_color=[[rowOddColor if i % 2 == 0 else rowEvenColor for i in range(df_plot.shape[0])]],
                    align=['left'] + ['center'] * (len(df_plot.columns)-1),
                    font=dict(color='darkslategray', size=11),
                    # format = ['html'] * len(df_plot.columns)  # Enable HTML formatting
                    )
    )
    ])
    table_fig.update_layout(
        annotations=[
            dict(
                text="*FLAG=TRUE: At least one kilobase-scale region overlapped by gene is <50% uniquely mappable or in the top 99.99th percentile of mutation rate.",
                x=0,
                y=-0.15,
                xref="paper",
                yref="paper",
                showarrow=False,
                align="left",
                valign="top",
                font=dict(size=12)
            )
        ]
    )

    return df_kept, pvals, pval_bounds, logfc, logq, logq_bounds, labels, ind_kept, table_fig

def generate_dig_non_coding_plots(df, mut_key, bur_key, display_bounds_key, scatterpoint_key, alp=0.1):
    # prepare plot data for all combinations of mut_type and burden_type dropdown options
    plot_data = {}
    mut_val = combined_mutation_type[mut_key]
    bur_val = non_coding_region_burden_type[bur_key]
    display_bounds_val = display_bounds_type[display_bounds_key] 
    scatterpoint_val = scatterpoint_type[scatterpoint_key]
    text_special = ""
   
    if not (mut_key in ['indels', 'indels_snvs', 'mutations'] and bur_key == 'sample_wise'):
        df_kept, pvals, pval_bounds, logfc, logq, logq_bounds, labels, ind_kept, table_fig = generate_non_coding_region_plot_data(df, mut_val, bur_val, display_bounds_val, scatterpoint_val)

        # Volcano Plot

        # Scatter plots
        ind_capped = logq > ymax
        ind_ncapped = np.logical_and(ind_kept, ~ind_capped)
        labels_capped = labels[ind_capped].tolist()
        logfc_capped = logfc[ind_capped].tolist()
        logq_capped = logq[ind_capped].tolist()
        if display_bounds_val:
            logq_upper = logq_bounds[:, 0] - logq
            logq_lower = logq - logq_bounds[:, 1]
            dict_erry_sig = dict(
                type='data',
                symmetric=False,
                array=np.round(logq_upper[ind_ncapped], 3).tolist(),
                arrayminus=np.round(logq_lower[ind_ncapped], 3).tolist(),
                thickness=thk_err,
                width=wid_err,
                color=col_err_sig
            )
            dict_erry_nonsig = dict(
                type='data',
                symmetric=False,
                array=np.round(logq_upper[~ind_kept], 3).tolist(),
                arrayminus=np.round(logq_lower[~ind_kept], 3).tolist(),
                thickness=thk_err,
                width=wid_err,
                color=col_err_nonsig
            )
            logq_upper_capped = logq_upper[ind_capped].tolist()
            logq_lower_capped = logq_lower[ind_capped].tolist()
            labels_capped = [(
                f"({logfc_capped[i]:.3f}, {logq_capped[i]:.3f} +{logq_upper_capped[i]:.3f} / -{logq_lower_capped[i]:.3f})<br>{l}")
                for i, l in enumerate(labels_capped)]
            ylim_upper = min(np.max(logq_upper + logq), ymax) * (1 + hor_buffer)
        else:
            dict_erry_sig, dict_erry_nonsig = None, None
            labels_capped = [(
                f"({logfc_capped[i]:.3f}, {logq_capped[i]:.3f})<br>{l}")
                for i, l in enumerate(labels_capped)]
            ylim_upper = min(np.max(logq), ymax) * (1 + hor_buffer)

        volcano_fig = go.Figure()
        volcano_fig.add_trace(
            go.Scatter(
                x=logfc[~ind_kept].tolist(),
                y=logq[~ind_kept].tolist(),
                error_y=dict_erry_nonsig,
                mode='markers',
                marker=dict(color=col_nonsig, opacity=opac_nonsig),
                text=labels[~ind_kept].tolist(),
                name='Non-significant',
                showlegend=False,
                xhoverformat='.3f',
                yhoverformat='.3f'
            )
        )

        volcano_fig.add_trace(
            go.Scatter(
                x=logfc[ind_ncapped].tolist(),
                y=logq[ind_ncapped].tolist(),
                error_y=dict_erry_sig,
                mode='markers',
                marker=dict(color=col_sig, opacity=opac_sig),
                text=labels[ind_ncapped].tolist(),
                name='Significant',
                showlegend=False,
                xhoverformat='.3f',
                yhoverformat='.3f'
            )
        )

        volcano_fig.add_trace(
            go.Scatter(
                x=logfc[ind_capped].tolist(),
                y=[ymax] * sum(ind_capped),
                mode='markers',
                marker=dict(color=col_sig, opacity=opac_sig),
                text=labels_capped,
                name='Significant',
                showlegend=False,
                hoverinfo='name+text',
            )
        )

        # Line plots
        volcano_fig.add_trace(
            go.Scatter(
                x=[1, 1],
                y=[0, ylim_upper],
                mode='lines',
                line=dict(dash=typ_thin, color=col_thin, width=thk_thin),
                showlegend=False,
                hoverinfo='skip'
            )
        )
        volcano_fig.add_trace(
            go.Scatter(
                x=[0, np.max(logfc) * (1 + hor_buffer)],
                y=[-np.log10(alp), -np.log10(alp)],
                mode='lines',
                line=dict(dash=typ_thick, color=col_thick, width=thk_thick),
                showlegend=False,
                hoverinfo='skip'
            )
        )
        if ymax < ylim_upper:
            volcano_fig.add_trace(
                go.Scatter(
                    x=[0, np.max(logfc) * (1 + hor_buffer)],
                    y=[ymax] * 2,
                    mode='lines',
                    line=dict(dash=typ_thin, color=col_thin, width=thk_thin),
                    showlegend=False,
                    hoverinfo='skip'
                )
            )

        # Formatting the figure
        volcano_fig.update_layout(
            title='Observed/Expected counts vs. False Discovery Rate:',
            xaxis_title='Log2(Observed/Expected + 1)',
            yaxis_title='-Log10(FDR)',
            xaxis=dict(range=[0, np.max(logfc) * (1 + hor_buffer)]),
            yaxis=dict(range=[0, ylim_upper]),
            template='plotly_white'
        )

        # Q-Q Plot

        # Scatter plots
        x = -np.log10(np.arange(1, len(pvals) + 1) / (len(pvals) + 1))
        y = -np.log10(pvals)
        ind_capped = y > ymax
        ind_ncapped = np.logical_and(ind_kept, ~ind_capped)
        labels_capped = labels[ind_capped].tolist()
        x_capped = x[ind_capped].tolist()
        y_capped = y[ind_capped].tolist()

        if display_bounds_val:
            y_upper = -np.log10(pval_bounds[:, 0]) - y
            y_lower = y + np.log10(pval_bounds[:, 1])
            dict_erry_sig = dict(
                type='data',
                symmetric=False,
                array=np.round(y_upper[ind_ncapped], 3).tolist(),
                arrayminus=np.round(y_lower[ind_ncapped], 3).tolist(),
                thickness=thk_err,
                width=wid_err,
                color=col_err_sig
            )
            dict_erry_nonsig = dict(
                type='data',
                symmetric=False,
                array=np.round(y_upper[~ind_kept], 3).tolist(),
                arrayminus=np.round(y_lower[~ind_kept], 3).tolist(),
                thickness=thk_err,
                width=wid_err,
                color=col_err_nonsig
            )
            y_upper_capped = y_upper[ind_capped].tolist()
            y_lower_capped = y_lower[ind_capped].tolist()
            labels_capped = [(
                f"({x_capped[i]:.3f}, {y_capped[i]:.3f} +{y_upper_capped[i]:.3f} / -{y_lower_capped[i]:.3f})<br>{l}")
                for i, l in enumerate(labels_capped)]
            ylim_upper = min(np.max(y_upper + y), ymax) * (1 + hor_buffer)
        else:
            dict_erry_sig, dict_erry_nonsig = None, None
            labels_capped = [(
                f"({x_capped[i]:.3f}, {y_capped[i]:.3f})<br>{l}")
                for i, l in enumerate(labels_capped)]
            ylim_upper = min(np.max(y), ymax) * (1 + hor_buffer)

        xi = np.arange(1, len(pvals) + 1)
        clower = -np.log10(sp.stats.beta.ppf((1 - ci) / 2, xi, xi[::-1]))
        cupper = -np.log10(sp.stats.beta.ppf((1 + ci) / 2, xi, xi[::-1]))

        qq_fig = go.Figure()
        qq_fig.add_trace(
            go.Scatter(
                x=x[~ind_kept].tolist(),
                y=y[~ind_kept].tolist(),
                error_y=dict_erry_nonsig,
                mode='markers',
                marker=dict(color=col_nonsig, opacity=opac_nonsig),
                text=labels[~ind_kept].tolist(),
                name='Non-significant',
                showlegend=False,
                xhoverformat='.3f',
                yhoverformat='.3f'
            )
        )

        qq_fig.add_trace(
            go.Scatter(
                x=x[ind_ncapped].tolist(),
                y=y[ind_ncapped].tolist(),
                error_y=dict_erry_sig,
                mode='markers',
                marker=dict(color=col_sig, opacity=opac_sig),
                text=labels[ind_ncapped].tolist(),
                name='Significant',
                showlegend=False,
                xhoverformat='.3f',
                yhoverformat='.3f'
            )
        )

        qq_fig.add_trace(
            go.Scatter(
                x=x_capped,
                y=[ymax] * sum(ind_capped),
                mode='markers',
                marker=dict(color=col_sig, opacity=opac_sig),
                text=labels_capped,
                name='Significant',
                showlegend=False,
                hoverinfo='name+text',
            )
        )

        # Line plots
        qq_fig.add_trace(
            go.Scatter(
                x=[0, np.max(x)],
                y=[0, np.max(x)],
                mode='lines',
                line=dict(dash=typ_thick, color=col_thick, width=thk_thick),
                showlegend=False,
                hoverinfo='skip'
            )
        )
        if ymax < ylim_upper:
            qq_fig.add_trace(
                go.Scatter(
                    x=[0, np.max(x) * (1 + hor_buffer)],
                    y=[ymax] * 2,
                    mode='lines',
                    line=dict(dash=typ_thin, color=col_thin, width=thk_thin),
                    showlegend=False,
                    hoverinfo='skip'
                )
            )

        # Area of confidence intervals for the identity line
        qq_fig.add_trace(go.Scatter(
            x=np.concatenate([x, x[::-1]]).tolist(),  # Combine x values for fill
            y=np.concatenate([clower, cupper[::-1]]).tolist(),  # Combine y values for fill
            fill='toself',
            fillcolor=col_beta,
            line=dict(color='rgba(255,255,255,0)'),  # No line for the filled area
            showlegend=False,
            hoverinfo='skip'
        ))

        # Formatting the figure
        qq_fig.update_layout(
            title='QQ-Plot of P-values:',
            xaxis_title='Expected -Log10(P-value)',
            yaxis_title='Observed -Log10(P-value)',
            xaxis=dict(range=[0, np.max(x) * (1 + hor_buffer)]),
            yaxis=dict(range=[0, ylim_upper]),
            template='plotly_white'
        )

        # Save figures as separate data
        plot_data[f"{mut_key}_{bur_key}_{display_bounds_key}_{scatterpoint_key}"] = {
            'volcano': volcano_fig.to_dict(),
            'qq': qq_fig.to_dict(),
            'table': table_fig.to_dict(),
            'text': bur_key + ' Mutation Burden of ' + mut_key,
            'textcolor': 'black-text'
        }
    else:
        plot_data[f"{mut_key}_{bur_key}_{display_bounds_key}_{scatterpoint_key}"] = {
            'volcano': None,
            'qq': None,
            'dnds': None,
            'table': None,
            'text': SPECIAL_TEXT,
            'textcolor': 'red-text'
        }
        text_special = SPECIAL_TEXT

    # # convert plot data to JSON-like structure
    # plot_data_json = json.dumps(plot_data)

    # generate static figures for the default values
    fig_mu = px.histogram(df_kept,
                          x='MU',
                          labels={'MU': 'MU'},
                          opacity=opac_bar,
                          log_y=True,
                          color_discrete_sequence=[col_bar])
    fig_mu.update_layout(
        title='Mean of GP model:',
        xaxis_title='MU (mutations per kilobase)',
        yaxis_title='Number of genes',
        template='plotly_white')

    fig_sigma = px.histogram(df_kept,
                             x='SIGMA',
                             labels={'SIGMA': 'SIGMA'},
                             opacity=opac_bar,
                             log_y=True,
                             color_discrete_sequence=[col_bar])
    fig_sigma.update_layout(
        title='Standard deviation of GP model:',
        xaxis_title='SIGMA (mutations per kilobase)',
        yaxis_title='Number of genes',
        template='plotly_white')

    return df_kept, volcano_fig, qq_fig, fig_mu, fig_sigma, table_fig, text_special

# dNdScv REPORT CODE

# dropdown options: number of genes to display in the scatter plots
# minimum number of rows (genes) to display in the table
n_rows_min = 50
# horizontal buffer for the scatter plots
hor_buffer = 0.05
hor_buffer_global = 0.25
# vertical gap between gene labels and scatter points
y_gap_annot = 0.01
# buffer for the number of rows in the table
n_rows_buffer = 0.5
# maximum value along the axis of p-values
pval_max = 16

# properties of significant points
col_sig = 'rgba(255, 0, 0, 1)'
opac_sig = 0.7
# properties of near-significant points
col_nearsig = 'rgba(30, 144, 255, 1)'
opac_nearsig = 0.7
# properties of non-significant points
col_nonsig = 'rgba(0, 0, 0, 1)'
opac_nonsig = 0.5
# properties of thinner lines
thk_thin = 0.75
col_thin = 'gray'
typ_thin = 'dash'
# properties of thicker lines
thk_thick = 2.5
col_thick = 'gray'
typ_thick = 'dash'
# properties of error bars
wid_err = 1
wid_err_global = 2
thk_err = 0.5
thk_err_global = 1
opac_err = 0.3
# font size of axes labels
fs_axes = 14

dnds_markers = {
    'mis': ['Missense', 'circle'],
    'trunc': ['Truncating', 'square'],
    'ind': ['InDel', 'star']
}

ngenes = {
    'All': 'all',
    'First 30': 30,
    'First 20': 20,
    'First 10': 10,
    'None': None
}

# derived parameters
col_err_sig = ','.join(col_sig.split(',')[:-1]) + ', {})'.format(opac_err)
col_err_nearsig = ','.join(col_nearsig.split(',')[:-1]) + ', {})'.format(opac_err)
col_err_nonsig = ','.join(col_nonsig.split(',')[:-1]) + ', {})'.format(opac_err)


def reformat_numbers(x, format='{:.2E}'):
    """
    Reformat numbers in an array to a specific format
    """
    return [format.format(n) for n in x]


def get_capped_labels(x_capped, y_capped, labels_capped):
    labels_capped = [(
        f"({x_capped[i]:.2f}, {y_capped[i]:.2f})<br>{l}")
        for i, l in enumerate(labels_capped)]
    return np.array(labels_capped)

def prepare_dnds_scatter_data(labels, x, y, y_lower, y_upper, ind, col_err, thk_err, wid_err, capped=False):
    """
    Prepare scatter plot data and error bars for given indices.

    Parameters:
    labels (np.ndarray): Array of labels.
    x (np.ndarray): Array of x values.
    y (np.ndarray): Array of y values.
    y_lower (np.ndarray): Array of lower y values.
    y_upper (np.ndarray): Array of upper y values.
    ind (np.ndarray): Boolean array indicating the indices to use.
    col_err (str): Color for the error bars.
    thk_err (float): Thickness of the error bars.
    wid_err (float): Width of the error bars.
    labels_selected (np.ndarray, optional): Array of selected labels. Defaults to None.

    Returns:
    dict: Dictionary containing scatter plot data and error bars.
    """
    x_selected = x[ind]
    y_selected = y[ind]
    y_selected_lower = y_lower[ind]
    y_selected_upper = y_upper[ind]
    y_selected_upper_err = y_selected_upper - y_selected
    y_selected_lower_err = y_selected - y_selected_lower
    if capped:
        labels_selected = labels[ind]
        labels_selected = np.array([(
            f"({x_selected[i]:.2f}, {y_selected[i]:.2f} +{y_selected_upper_err[i]:.2f} / -{y_selected_lower_err[i]:.2f})<br>{l}")
            for i, l in enumerate(labels_selected)])
        x[ind] = pval_max
        x_selected = np.array([pval_max] * ind.sum())
    else:
        labels_selected = labels[ind]
    dict_erry = dict(
        type='data',
        symmetric=False,
        array=np.round(y_selected_upper_err, 2).tolist(),
        arrayminus=np.round(y_selected_lower_err, 2).tolist(),
        thickness=thk_err,
        width=wid_err,
        color=col_err
    )
    return {
        'labels': labels_selected,
        'x': x_selected,
        'y': y_selected,
        'dict_erry': dict_erry
    }


def plot_dnds(df, mut_typ, mut_name, alp, alp_nearsig):
    """
    Generate a scatter plot with FDR values on the horizontal and dNdS values on the vertical axis.

    Parameters:
    df (pd.DataFrame): DataFrame containing dNdS values and confidence intervals.
    mut_typ (str): Mutation type (e.g., 'mis' for missense, 'tru' for truncating).
    mut_name (str): Mutation name for labeling the plot.
    alp (float): Significance level for determining significant genes.
    alp_nearsig (float): Near-significance threshold.

    Returns:
    plotly.graph_objects.Figure: Plotly figure object containing the scatter plot.
    """
    labels = df.index.copy().to_numpy()
    x = -np.log10(df.qglobal_cv.copy().to_numpy())
    y = df[mut_typ + '_mle'].copy().to_numpy()
    y_lower = df[mut_typ + '_low'].copy().to_numpy()
    y_upper = df[mut_typ + '_high'].copy().to_numpy()
    # indicator of significant genes
    ind_sig = x >= -np.log10(alp)
    # indicator of capped FDR values
    ind_capped = x > pval_max
    # indicator of significant genes that are not capped
    ind_ncapped = np.logical_and(ind_sig, ~ind_capped)
    # indicator of significant or near-significant genes
    ind_hits = x >= -np.log10(alp_nearsig)
    # indicator of near-significant genes
    ind_nearsig = np.logical_and(~ind_sig, ind_hits)
    # indicator of non-significant genes
    ind_nsig = ~ind_sig

    # scatter points for non-significant genes
    nsig_data = prepare_dnds_scatter_data(labels, x, y, y_lower, y_upper, ind_nsig, col_err_nonsig, thk_err, wid_err)
    # scatter points for nearly-significant genes
    nearsig_data = prepare_dnds_scatter_data(labels, x, y, y_lower, y_upper, ind_nearsig, col_err_nearsig, thk_err, wid_err)
    # scatter points for significant genes with non-capped FDR values
    ncapped_data = prepare_dnds_scatter_data(labels, x, y, y_lower, y_upper, ind_ncapped, col_err_sig, thk_err, wid_err)
    # scatter points for significant genes with capped FDR values
    capped_data = prepare_dnds_scatter_data(labels, x, y, y_lower, y_upper, ind_capped, col_err_sig, thk_err, wid_err, capped=True)

    plotly_fig = go.Figure()

    scatter_data = [
        (nsig_data['x'].tolist(), nsig_data['y'].tolist(), nsig_data['dict_erry'], col_nonsig, opac_nonsig,
         nsig_data['labels'].tolist(), 'Non-significant'),
        (nearsig_data['x'].tolist(), nearsig_data['y'].tolist(), nearsig_data['dict_erry'], col_nearsig, opac_nearsig,
         nearsig_data['labels'].tolist(), 'Nearly-significant'),
        (ncapped_data['x'].tolist(), ncapped_data['y'].tolist(), ncapped_data['dict_erry'], col_sig, opac_sig,
         ncapped_data['labels'].tolist(), 'Significant'),
        (capped_data['x'].tolist(), capped_data['y'].tolist(), capped_data['dict_erry'], col_sig, opac_sig,
         capped_data['labels'], 'Significant')
    ]

    for i, (x_data, y_data, error_y, color, opacity, text, name) in enumerate(scatter_data):
        if i == len(scatter_data)-1:
            hoverinfo = 'name+text'
        else:
            hoverinfo = 'all'
        plotly_fig.add_trace(
            go.Scatter(
                x=x_data,
                y=y_data,
                error_y=error_y,
                mode='markers',
                marker=dict(color=color, opacity=opacity),
                text=text,
                name=name,
                showlegend=False,
                hoverinfo=hoverinfo,
                xhoverformat='.2f',
                yhoverformat='.2f'
            )
        )

    # define boundaries
    xlim = [0, max(x) * (1 + hor_buffer)]
    ylim = [0, max(max(y_upper[ind_sig]), max(y)) * (1 + hor_buffer)]
    # plot cap line
    plotly_fig.add_trace(
        go.Scatter(
            x=[pval_max] * 2,
            y=ylim,
            mode='lines',
            line=dict(dash=typ_thin, color=col_thin, width=thk_thin),
            showlegend=False,
            hoverinfo='skip'
        )
    )
    # format and add lables
    plotly_fig.update_layout(
        title=dict(
            text='dNdS' + ' Ratio of ' + (mut_name + ' mutations:').title(),
            font=dict(size=20)  # title font size
        ),
        xaxis_title='-Log10(FDR)',
        yaxis_title='dNdS ratio',
        xaxis=dict(
            range=xlim
        ),
        yaxis=dict(
            range=ylim
        ),
        template='plotly_white'
    )
    return plotly_fig


def plot_dnds_qq(x, y, labels, col, opac, mark, name, fig, hi='x+y+text+name'):
    
    fig.add_trace(
        go.Scatter(
            x=x.tolist(),
            y=y.tolist(),
            mode='markers',
            marker=dict(color=col, opacity=opac, symbol=mark, size=msize),
            text=labels.tolist(),
            name=name,
            hoverinfo=hi,
            showlegend=False,
            xhoverformat='.2f',
            yhoverformat='.2f'
        )
    )
    return fig


def generate_dnds_dataframe(
        path_dnds_out,
        path_dnds_ci,
        path_dnds_global,
):
    """
    Generate a dNdScv report for coding regions.

    Parameters:
    path_dnds_out (str): Path to the dNdScv output file containing p-values.
    path_dnds_ci (str): Path to the dNdScv confidence intervals file.
    path_dnds_global (str): Path to the dNdScv global values file.
    path_cgc_list (str): Path to the list of CGC genes.
    path_pancan_list (str): Path to the list of PanCanAtlas genes.
    dir_output (str): Output directory.
    prefix_output (str, optional): Prefix for the output file. Defaults to None.
    alp (float, optional): Significance level. Defaults to 0.1.
    alp_nearsig (float, optional): Near-significance threshold (default: 0.25).
    Returns:
    None
    """
    # Driver gene lists
    cgc_list_path = "gs://getzlab-workflows-reference_files-oa/hg19/dig/cancer_gene_census_2024_06_20.tsv"
    pancan_list_path = "gs://getzlab-workflows-reference_files-oa/hg19/dig/pancanatlas_genes.tsv"

    # lists of known driver genes
    cgc_list = pd.read_csv(cgc_list_path, sep='\t').to_numpy().flatten()
    pancan_list = pd.read_csv(pancan_list_path, sep='\t').to_numpy().flatten()
    # output from dNdScv containing p-values
    df_out = pd.read_csv(path_dnds_out, sep='\t').set_index('gene_name')
    # output from dNdScv containing dNdS values and confidence intervals
    df_ci = pd.read_csv(path_dnds_ci, sep='\t').set_index('gene')
    # merging the two dataframes
    df_merged = pd.concat([df_out, df_ci], axis=1)
    # global dNdS values
    df_global = pd.read_csv(path_dnds_global, sep='\t')
    # dataframe that will be plotted in a table format
    df_plot = df_merged.copy().reset_index().rename(columns={'index': 'GENE'})
    df_plot['RANK'] = df_plot.index + 1
    # adding indicator of genes being part of the CGC or PanCan list
    df_plot['CGC'] = df_plot.GENE.isin(cgc_list)
    df_plot['PANCAN'] = df_plot.GENE.isin(pancan_list)

    # renaming columns
    dict_n = {c: c.upper() for c in ['n_syn', 'n_mis', 'n_non', 'n_spl', 'n_ind']}
    dict_w = {c: 'dNdS_' + c.split('_')[0][1:].upper() for c in ['wmis_cv', 'wnon_cv', 'wspl_cv', 'wind_cv']}
    dict_pvals = {c: 'PVAL_' + c.split('_')[0][1:].upper() for c in ['pmis_cv', 'ptrunc_cv', 'pind_cv']}
    df_plot = df_plot.rename(columns=dict_n)
    df_plot = df_plot.rename(columns=dict_w)
    df_plot = df_plot.rename(columns=dict_pvals)
    df_plot = df_plot.rename(columns={'pglobal_cv': 'PVAL', 'qglobal_cv': 'FDR'})

    return df_plot, df_merged, df_global

def generate_dnds_report(
        df_plot, 
        df_merged,
        df_global,
        ngenes_key,
        alp=0.1,
        alp_nearsig=0.25):
    """ 
    """
    dict_n = {c: c.upper() for c in ['n_syn', 'n_mis', 'n_non', 'n_spl', 'n_ind']}
    dict_w = {c: 'dNdS_' + c.split('_')[0][1:].upper() for c in ['wmis_cv', 'wnon_cv', 'wspl_cv', 'wind_cv']}
    dict_pvals = {c: 'PVAL_' + c.split('_')[0][1:].upper() for c in ['pmis_cv', 'ptrunc_cv', 'pind_cv']}

    cols_kept = (['RANK', 'GENE'] + list(dict_n.values()) + list(dict_w.values()) + list(dict_pvals.values()) +
                 ['PVAL', 'FDR', 'CGC', 'PANCAN'])
    # cols_kept = (['RANK', 'GENE'] + list(dict_n.values()) + list(dict_w.values()) + ['PVAL', 'FDR', 'CGC', 'PANCAN'])

    # # generate the dropdown options
    # ngenes_options = "\n".join([f'<option value="{key}">{key}</option>' for key in ngenes.keys()])


    # 'coding': ["Coding", "circle"],
    # 'promoter': ["Promoter", "square"],
    # '5utr': ["5\' UTR", "star"],
    # '3utr': ["3\' UTR", "triangle-down"]

    tests_pvals = ['p' + mut_typ + '_cv' for mut_typ in dnds_markers.keys()]
    # tests_pvals = [mut_typ + '_cv' for mut_typ in markers.keys()]
    # mutations = ['mis',	'trunc', 'allsubs',	'ind'] # looked at the dnds_out.tsv file to find the mutation types
    # tests_pvals = ["p" + mut_type + "_cv" for mut_type in mutations]

    # prepare plot data for all dropdown options
    plot_data = {}
    print("checking ngenes_key: ", ngenes_key)
    ngenes_val = ngenes[ngenes_key]
    # for ngenes_key, ngenes_val in ngenes.items():

    # Q-Q Plot
    # scatter plots
    labels = df_merged.index.copy().to_numpy()
    pvals = df_merged.pglobal_cv.copy().to_numpy()
    x = -np.log10(np.arange(1, len(pvals) + 1) / (len(pvals) + 1))
    y = -np.log10(pvals)
    logfdr = -np.log10(df_merged.qglobal_cv.copy().to_numpy())
    # indicator of "dominant" test (with minimal FDR) for each gene
    # test_dom = df_plot[tests_pvals].idxmin(axis=1).to_numpy()
    print("test_pvals: ", tests_pvals)
    test_dom = df_merged[tests_pvals].idxmin(axis=1).to_numpy()

    # indicator of genes whose p-value exceeds the maximum value
    ind_capped = y > pval_max
    # indicator of significant genes
    ind_sig = logfdr >= -np.log10(alp)
    # indicator of significant or near-significant genes
    ind_hits = logfdr >= -np.log10(alp_nearsig)
    # indicator of near-significant genes
    ind_nearsig = np.logical_and(~ind_sig, ind_hits)
    # indicator of non-significant genes
    ind_nsig = ~ind_sig
    # indicator of hits that are not capped
    ind_ncapped = np.logical_and(ind_sig, ~ind_capped)
    # indicator of hits that are not capped
    ind_ncapped_hits = np.logical_and(ind_hits, ~ind_capped)
    ylim_upper = min(np.max(y), pval_max) * (1 + hor_buffer)

    # Scatter plot
    qq_fig = go.Figure()

    for i in range(len(tests_pvals)):
        test_i = tests_pvals[i][1:-3]
        ind = test_dom == tests_pvals[i]
        i_nsig = np.logical_and(ind, ind_nsig)
        i_nearsig = np.logical_and(ind, ind_nearsig)
        i_ncapped = np.logical_and(ind, ind_ncapped)
        i_capped = np.logical_and(ind, ind_capped)

        qq_fig = plot_dnds_qq(x[i_nsig], y[i_nsig], labels[i_nsig], col_nonsig, opac_nonsig, dnds_markers[test_i][1],
                            'Non-significant', qq_fig)
        qq_fig = plot_dnds_qq(x[i_nearsig], y[i_nearsig], labels[i_nearsig], col_nearsig, opac_nearsig,
                            dnds_markers[test_i][1],
                            'Near-significance', qq_fig)
        qq_fig = plot_dnds_qq(x[i_ncapped], y[i_ncapped], labels[i_ncapped], col_sig, opac_sig, dnds_markers[test_i][1],
                            'Significant', qq_fig)
        xi_capped = x[i_capped]
        qq_fig = plot_dnds_qq(xi_capped, np.array([pval_max] * sum(i_capped)),
                            get_capped_labels(xi_capped, y[i_capped], labels[i_capped]), col_sig, opac_sig,
                            dnds_markers[test_i][1], 'Significant', qq_fig, hi='name+text')

    # Dummy points for legend
    tests_dom, test_counts = np.unique(test_dom, return_counts=True)
    for test_pval_i in tests_dom[np.argsort(test_counts)[::-1]]:
        test_i = test_pval_i[1:-3]
        qq_fig.add_trace(
            go.Scatter(
                y=[None],
                mode='markers',
                marker=dict(
                    color='white',
                    symbol=dnds_markers[test_i][1],
                    size=msize * 1.25,
                    line=dict(color='black', width=2)
                ),
                name=dnds_markers[test_i][0]
            )
        )

    # Add rotated text labels as annotations
    x_capped = x[ind_capped]
    y_capped = y[ind_capped]
    labels_capped = get_capped_labels(x_capped, y_capped, labels[ind_capped])

    # if ngenes_val is not None:
    if ngenes_val != 'none':
        print("ngenes_val: ", ngenes_val)
        if ngenes_val == 'all':
            x_annot = [x_capped, x[ind_ncapped_hits]]
            y_annot = [[pval_max] * sum(ind_capped), y[ind_ncapped_hits]]
            labels_annot = [labels_capped, labels[ind_ncapped_hits]]
        else:
            if ngenes_val <= len(x_capped):
                x_annot = [x_capped[:ngenes_val], []]
                y_annot = [[pval_max] * ngenes_val, []]
                labels_annot = [labels_capped[:ngenes_val], []]
            else:
                end_ncapped = ngenes_val - len(x_capped)
                x_annot = [x_capped, x[ind_ncapped_hits][:end_ncapped]]
                y_annot = [[pval_max] * sum(ind_capped), y[ind_ncapped_hits][:end_ncapped]]
                labels_annot = [labels_capped, labels[ind_ncapped_hits][:end_ncapped]]
        count = 0
        for (xa, ya, laba) in zip(x_annot, y_annot, labels_annot):
            for (xi, yi, label) in zip(xa, ya, laba):
                if ind_nearsig[count]:
                    coli = col_nearsig
                else:
                    coli = col_sig
                qq_fig.add_annotation(
                    x=xi,
                    y=yi - ylim_upper * y_gap_annot,
                    text=label.split('<br>')[-1],
                    showarrow=False,
                    font=dict(color=coli),
                    textangle=-90,
                    xanchor="center",
                    yanchor="top"
                )
                count += 1

    # line plots
    qq_fig.add_trace(
        go.Scatter(
            x=[0, np.max(x) * (1 + hor_buffer)],
            y=[0, np.max(x) * (1 + hor_buffer)],
            mode='lines',
            line=dict(dash=typ_thick, color=col_thick, width=thk_thick),
            showlegend=False,
            hoverinfo='skip'
        )
    )
    qq_fig.add_trace(
        go.Scatter(
            x=[0, np.max(x) * (1 + hor_buffer)],
            y=[pval_max] * 2,
            mode='lines',
            line=dict(dash=typ_thin, color=col_thin, width=thk_thin),
            showlegend=False,
            hoverinfo='skip'
        )
    )
    # figure formatting
    qq_fig.update_layout(
        title='QQ-Plot of P-values:',
        xaxis_title='Expected -Log10(P-value)',
        yaxis_title='Observed -Log10(P-value)',
        xaxis=dict(range=[0, np.max(x) * (1 + hor_buffer)]),
        yaxis=dict(range=[0, ylim_upper]),
        template='plotly_white',
        legend=dict(
            title=dict(
                text="Dominant Mutation Type:",
            ),
            indentation=10
        )
    )
    # save figures as separate data
    plot_data[f"{ngenes_key}"] = {
        'qq': qq_fig.to_dict()
    }

    # Global dNdS Plot
    mle = df_global.mle
    upper_global = df_global.cihigh - mle
    lower_global = mle - df_global.cilow

    dict_erry_global = dict(
        type='data',
        symmetric=False,
        array=np.round(upper_global, 2).tolist(),
        arrayminus=np.round(lower_global, 2).tolist(),
        thickness=thk_err_global,
        width=wid_err_global,
        color=col_err_nonsig
    )

    x_labels = ['Missense', 'Nonsense', 'Splice site', 'Truncating', 'All']
    x_global = [i + 1 for i in range(len(df_global))]
    xlim_global = [-hor_buffer_global + min(x_global), max(x_global) + hor_buffer_global]

    fig_dnds_global = go.Figure()
    fig_dnds_global.add_trace(
        go.Scatter(
            x=xlim_global,
            y=[1, 1],
            mode='lines',
            line=dict(dash=typ_thick, color=col_thick, width=thk_thick),
            showlegend=False,
            hoverinfo='skip'
        )
    )
    fig_dnds_global.add_trace(
        go.Scatter(
            x=x_global,
            y=df_global['mle'].tolist(),
            error_y=dict_erry_global,
            mode='markers',
            marker=dict(color=col_nonsig, opacity=opac_nonsig),
            showlegend=False,
            xhoverformat='.2f',
            yhoverformat='.2f'
        )
    )

    fig_dnds_global.update_layout(
        title=dict(
            text='dNdS' + ' ratio across all mutations:'.title(),
            # font=dict(size=20)  # Set your desired title font size
        ),
        xaxis_title='',
        yaxis_title='dNdS ratio',
        xaxis=dict(
            range=xlim_global,
            tickvals=x_global,
            ticktext=x_labels,
            tickfont=dict(size=fs_axes)  # size of the tick labels
        ),
        yaxis=dict(
            titlefont=dict(size=fs_axes)  # size of y-axis title
        ),
        template='plotly_white'
    )

    #
    # Missense dNdS Plot
    fig_dnds_mis = plot_dnds(df_merged, 'mis', 'missense', alp, alp_nearsig)

    #
    # Truncating dNdS Plot
    fig_dnds_tru = plot_dnds(df_merged, 'tru', 'truncating', alp, alp_nearsig)

    #
    # Table Plot
    n_rows = max(n_rows_min, int(np.sum(ind_sig) * (1 + n_rows_buffer)))
    df_plot = df_plot.iloc[:n_rows][cols_kept].copy()

    for col in ['PVAL', 'FDR'] + list(dict_pvals.values()):
    # for col in ['PVAL', 'FDR']:
        df_plot[col] = reformat_numbers(df_plot[col].to_numpy())
    for col in list(dict_w.values()):
        df_plot[col] = reformat_numbers(df_plot[col].to_numpy(), format='{:.2f}')
    for c in list(dict_n.values()) + ['RANK']:
        df_plot[c] = df_plot[c].astype(int)

    headerColor = 'grey'
    rowEvenColor = 'lightgrey'
    rowOddColor = 'white'

    # making the significant rows bold
    df_plot = df_plot.astype(str)
    for i in range(df_plot.shape[0]):
        if ind_sig[i]:
            df_plot.loc[i, :] = '<b>' + df_plot.loc[i, :].astype(str) + '</b>'
            
    # generate table figure
    fig_table = go.Figure(data=[go.Table(
        header=dict(values=['<b>' + col + '</b>' for col in df_plot.columns],
                    line_color='darkslategray',
                    fill_color=headerColor,
                    align=['left'] + ['center'] * (len(df_plot.columns) - 1),
                    font=dict(color='white', size=12)
                    ),
        cells=dict(values=[df_plot[col].tolist() for col in df_plot.columns],
                   line_color='darkslategray',
                   fill_color=[[rowOddColor if i % 2 == 0 else rowEvenColor for i in range(df_plot.shape[0])]],
                   align=['left'] + ['center'] * (len(df_plot.columns) - 1),
                   font=dict(color='darkslategray', size=11),
                   )
                )
            ])
    
    return qq_fig, fig_dnds_global, fig_dnds_mis, fig_dnds_tru, df_plot
    # NEED TO RETURN ALL THE dNdS plots!!!



# CODE FOR COMPARISON dNdScv report
import argparse
import json
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from statsmodels.stats.multitest import fdrcorrection
import scipy as sp
from itertools import chain, combinations

#
# PARAMETERS

# Scatter Plot

# figure size
fig_size = (800, 700)
# fontsize
fs = 14
# padding for plot range
xy_pad = 0.05
# tick step size
tick_step = 2
# vertical gap between the annotation labels
ay_gap = 15
# colors
color_both = 'green'
color_only1 = 'magenta'
color_only2 = 'red'
color_neither = 'black'
# opacity
opacity=0.6

# Table Plot

headerColor = 'grey'
rowEvenColor = 'lightgrey'
rowOddColor = 'white'
lineColor = 'darkslategray'

# Static Table Plot
# minimum number of rows (genes) to display in the table
n_rows_min = 50
# buffer for the number of rows in the table
n_rows_buffer = 0.5

# minimum threshold at which we cap all FDRs
FDR_limit = 1e-16

log10FDR_limit = -np.log10(FDR_limit)

# dropdown options: methods to display results for
methods = {
    'MutSig2CV vs dNdScv': ['MutSig2CV', 'dNdScv'],
    'MutSig2CV vs DIG': ['MutSig2CV', 'DIG'],
    'dNdScv vs DIG': ['dNdScv', 'DIG']
}
#
# UTILITY FUNCTIONS

def reformat_numbers(x, format='{:.2E}'):
    """
    Reformat numbers in an array to a specific format
    """
    return [format.format(n) for n in x]

def format_cols(df):
    """
    Format columns in a DataFrame for better readability.

    This function reformats numerical columns containing 'PVAL' or 'FDR' in their names,
    converts 'RANK' and 'SIZE_coding' columns to integer strings, and replaces NaN values
    with 'NA' in all columns.

    Parameters
    ----------
    df : pandas.DataFrame
        The input DataFrame to be formatted.

    Returns
    -------
    pandas.DataFrame
        The formatted DataFrame with numerical columns reformatted, integer columns converted
        to strings, and NaN values replaced with 'NA'.
    """
    # formatting data in columns
    for col in [c for c in df.columns if 'PVAL' in c or 'FDR' in c]:
        df[col] = reformat_numbers(df[col].to_numpy())
    for c in ['RANK', 'CHROM', 'SIZE_coding']:
        is_nan = df[c].isna()
        v = np.zeros(df.shape[0], dtype='object')
        v[~is_nan] = df[c][~is_nan].astype(int).astype(str)
        v[is_nan] = 'NA'
        df[c] = v
    df = df.astype(str)
    for c in df.columns:
        df.loc[df[c].str.lower().isin(['nan', 'na']), c] = 'NA'
    return df

def find_overlapping_points(X, Y):
    """
    Identify overlapping points (points with identical x and y coordinates) in a dataset.

    Parameters:
    -----------
    X : array-like
        A 1D array of x-coordinates of the points.
    Y : array-like
        A 1D array of y-coordinates of the points.

    Returns:
    --------
    pts_olap : numpy.ndarray
        A 1D array of unique indices corresponding to points that overlap (i.e., have matching x and y coordinates).
    eq : list of numpy.ndarray
        A list of arrays, where each array contains the indices of a group of points that overlap.
        Each group represents a unique set of overlapping points.

    Notes:
    ------
    - The function compares every pair of points to check if their x and y coordinates match.
    - Groups of overlapping points are determined and stored such that no group is a subset of another.
    - If there are no overlapping points, `eq` will be an empty list, and `pts_olap` will be an empty array.

    Example:
    --------
    >>> X = [1, 2, 3, 1]
    >>> Y = [4, 5, 6, 4]
    >>> pts_olap, eq = find_overlapping_points(X, Y)
    >>> print(pts_olap)
    [0, 3]
    >>> print(eq)
    [array([0, 3])]
    """
    # matrix of booleans indicating whether two points (both their x and y coordinates) match
    ind_eq = np.logical_and(np.reshape(X, (-1, 1)) - X == 0, np.reshape(Y, (-1, 1)) - Y == 0)
    # collect all groups of matching points
    eq = []
    for i in range(len(X)):
        # indices of group of matching points
        eqi = set((np.where(ind_eq[i, i:])[0] + i).tolist())
        if len(eqi) > 1:  # check if there are at least two matching points
            # check whether the group of points are a subset of a previously stored larger group; if not, keep it
            addit = True
            for eqj in eq:
                if len(eqi - eqj) == 0:
                    addit = False
                    break
            if addit:
                eq.append(eqi)
    # indices of all points that overlap
    pts_olap = np.unique([item for sublist in eq for item in sublist])
    return pts_olap, [np.array(list(e)) for e in eq]

def preprocess_results(
        path_mutsig,
        path_dndscv,
        path_dig,
        # path_cgc,
        # path_pancan,
        alp=0.1
):
    """
    Preprocess and combine results from multiple statistical methods to identify significant genes.

    This function loads and formats output files from MutSig2CV, dNdScv, and DIG methods, combines their results,
    and identifies significant genes based on a specified significance threshold (alpha). It also adds indicators
    for Cancer Gene Census (CGC) and PanCanAtlas membership.

    Parameters
    ----------
    path_mutsig : str
        Path to the MutSig2CV output file (sig_genes.txt).
    path_dndscv : str
        Path to the dNdScv output file (dnds_out.tsv).
    path_dig : str
        Path to the DIG output file (combined.dig.results.txt).
    path_cgc : str
        Path to the Cancer Gene Census (CGC) reference file (cancer_gene_census.tsv).
    path_pancan : str
        Path to the PanCanAtlas genes reference file (pancanatlas_genes.tsv).
    alp : float, optional
        Significance threshold for the FDR values (default is 0.1).

    Returns
    -------
    pandas.DataFrame
        A DataFrame containing combined and processed results with the following columns:
        - 'RANK': Rank of the gene based on combined p-value.
        - 'GENE': Gene identifier.
        - 'SIZE_coding': Coding region size.
        - 'FDR_MutSig2CV': FDR from MutSig2CV.
        - 'FDR_dNdScv': FDR from dNdScv.
        - 'FDR_DIG': FDR from DIG.
        - 'FDR_min': Minimum value of FDR out of the three tests.
        - 'SIG_MutSig2CV': Significance indicator for MutSig2CV.
        - 'SIG_dNdScv': Significance indicator for dNdScv.
        - 'SIG_DIG': Significance indicator for DIG.
        - 'CGC': Indicator for CGC membership.
        - 'PANCAN': Indicator for PanCanAtlas membership.
    """
    # load MutSig2CV output
    df_mutsig = pd.read_csv(path_mutsig, sep='\t', low_memory=False).set_index('gene')
    df_mutsig = df_mutsig.rename(columns={'q': 'FDR_MutSig2CV', 'p': 'PVAL_MutSig2CV', 'codelen': 'SIZE_MutSig2CV'})[
        ['SIZE_MutSig2CV', 'PVAL_MutSig2CV', 'FDR_MutSig2CV']]
    # load dNdScv output
    df_dndscv = pd.read_csv(path_dndscv, sep='\t').set_index('gene_name')
    df_dndscv = df_dndscv.rename(columns={'qglobal_cv': 'FDR_dNdScv', 'pglobal_cv': 'PVAL_dNdScv'})[
        ['PVAL_dNdScv', 'FDR_dNdScv']]
    # load DIG output
    df_dig = pd.read_csv(path_dig, sep='\t', low_memory=False).set_index('GENE')
    df_dig = df_dig.rename(columns={'SIZE_coding': 'SIZE_DIG', 'PVAL_coding_MUT_recalc': 'PVAL_DIG'})[
        ['CHROM', 'SIZE_DIG', 'PVAL_DIG']]
    df_dig['FDR_DIG'] = fdrcorrection(df_dig['PVAL_DIG'].astype(float))[1]
    df_dig['SIZE_DIG'] = df_dig['SIZE_DIG'].astype(int)
    # add coding region sizes
    df = pd.concat([df_mutsig, df_dndscv, df_dig], axis=1)
    df['SIZE_coding'] = df['SIZE_MutSig2CV'].copy()
    is_nan = df['SIZE_coding'].isna().copy()
    df.loc[is_nan, 'SIZE_coding'] = df['SIZE_DIG'][is_nan].copy()
    # combine p-values
    cols_pval = [c for c in df.columns if 'PVAL' in c]
    for col in cols_pval:
        method = col.split('_')[-1]
        is_nan = df[f'FDR_{method}'].isna()
        df.loc[~is_nan, f'SIG_{method}'] = df.loc[~is_nan, f'FDR_{method}'] < alp
    df['FDR_min'] = df[['FDR_MutSig2CV', 'FDR_dNdScv', 'FDR_DIG']].min(axis=1)
    df = df.sort_values('FDR_min').reset_index().rename(columns={'index': 'GENE'})
    df['RANK'] = df.index + 1
    df = df[['RANK', 'GENE', 'CHROM', 'SIZE_coding', 'FDR_MutSig2CV', 'FDR_dNdScv', 'FDR_DIG', 'FDR_min', 'SIG_MutSig2CV', 'SIG_dNdScv', 'SIG_DIG']]
    

    path_cgc = "gs://getzlab-workflows-reference_files-oa/hg19/dig/cancer_gene_census_2024_06_20.tsv"
    path_pancan = "gs://getzlab-workflows-reference_files-oa/hg19/dig/pancanatlas_genes.tsv"
    # path_cgc = ''
    # path_pancan = ''
    # add indicators of CGC and PanCanAtlas membership
    cgc_list = pd.read_csv(path_cgc, sep='\t').to_numpy().flatten()
    pancan_list = pd.read_csv(path_pancan, sep='\t').to_numpy().flatten()
    df['CGC'] = df.GENE.isin(cgc_list).copy()
    df['PANCAN'] = df.GENE.isin(pancan_list).copy()

    return df
#
# PLOTLY FIGURE/TABLE GENERATING FUNCTIONS

def plot_fdr_comparison(df, method1, method2, alp):
    """
    Generate a scatter plot using Plotly to compare -log10(FDR) values between two statistical methods.

    This function categorizes genes into four significance groups:
    1. Significant in both methods.
    2. Significant in only the first method.
    3. Significant in only the second method.
    4. Not significant in either method.

    Overlapping points are identified and annotated, FDR values are capped at a predefined threshold,
    and plot features such as identity lines and significance thresholds are highlighted.

    Parameters
    ----------
    df : pandas.DataFrame
        Input DataFrame containing FDR values for the two methods. The DataFrame must include:
        - 'FDR_{method1}' and 'FDR_{method2}': Columns with FDR values for each method.
        - 'GENE': Column with gene identifiers.

    method1 : str
        Name of the first method used for FDR calculation. This is used for column references and axis labels.

    method2 : str
        Name of the second method used for FDR calculation. This is used for column references and axis labels.

    alp : float
        Significance threshold for the FDR values. Genes with FDR values below this threshold are considered significant.

    Returns
    -------
    tuple
        A tuple containing:
        1. fig : plotly.graph_objects.Figure
            The generated scatter plot figure.
        2. sig_genes_dict : dict
            A dictionary with keys representing three significance categories:
            - `'{method1} and {method2}'`: Genes significant in both methods.
            - `'{method1} only'`: Genes significant only in the first method.
            - `'{method2} only'`: Genes significant only in the second method.
            Each key contains a NumPy array of gene identifiers.

    Notes
    -----
    - The FDR values are transformed to -log10(FDR) for visualization.
    - FDR values exceeding a predefined threshold are capped for scale consistency.
    - Overlapping points are grouped and annotated to minimize plot clutter.

    Example
    -------
    >>> result = plot_fdr_comparison(df, method1='Method_A', method2='Method_B', alp=0.1)
    >>> print(result[1]['Method_A and Method_B'])
    ['Gene1', 'Gene2', 'Gene3']
    """

    log10alp = -np.log10(alp)

    #
    # ASSEMBLE DATA FRAME

    # assemble data frame that will be plotted
    df_plot = df[['GENE']].copy()
    # original values
    df_plot[f"{method1}: -log10(FDR)"] = -np.log10(df[f'FDR_{method1}'])
    df_plot[f"{method2}: -log10(FDR)"] = -np.log10(df[f'FDR_{method2}'])
    # capped values
    df_plot[f"{method1}: -log10(FDR) capped"] = df_plot[f"{method1}: -log10(FDR)"].copy()
    df_plot[f"{method2}: -log10(FDR) capped"] = df_plot[f"{method2}: -log10(FDR)"].copy()
    df_plot.loc[
        df_plot[f"{method1}: -log10(FDR) capped"] > log10FDR_limit, f"{method1}: -log10(FDR) capped"] = log10FDR_limit
    df_plot.loc[
        df_plot[f"{method2}: -log10(FDR) capped"] > log10FDR_limit, f"{method2}: -log10(FDR) capped"] = log10FDR_limit
    # annotate genes based on their significance with either of the two methods
    df_plot.loc[np.logical_and(df_plot[f"{method1}: -log10(FDR)"] >= log10alp,
                               df_plot[f"{method2}: -log10(FDR)"] >= log10alp), 'Significant'] = 'both'
    df_plot.loc[np.logical_and(df_plot[f"{method1}: -log10(FDR)"] >= log10alp,
                               df_plot[f"{method2}: -log10(FDR)"] < log10alp), 'Significant'] = f'{method1} only'
    df_plot.loc[np.logical_and(df_plot[f"{method1}: -log10(FDR)"] < log10alp,
                               df_plot[f"{method2}: -log10(FDR)"] >= log10alp), 'Significant'] = f'{method2} only'
    df_plot.loc[np.logical_and(df_plot[f"{method1}: -log10(FDR)"] < log10alp,
                               df_plot[f"{method2}: -log10(FDR)"] < log10alp), 'Significant'] = 'neither'

    #
    # PRE-PROCESS DATA

    # displayed range
    xy_max = max([df_plot[f"{method1}: -log10(FDR) capped"].max(), df_plot[f"{method2}: -log10(FDR) capped"].max()]) * (
                1 + xy_pad)

    #
    # ASSEMBLE PLOTLY FIGURE

    fig = go.Figure()
    # plot "identity" line
    fig.add_trace(go.Scatter(x=[0, xy_max], y=[0, xy_max], showlegend=False, hoverinfo='skip', mode='lines',
                             line=dict(color='gray', dash='dash')))
    # plot lines that indicate significance thresholds
    fig.add_trace(go.Scatter(x=[0, xy_max], y=[log10alp] * 2, showlegend=False, hoverinfo='skip', mode='lines',
                             line=dict(color='gray', dash='solid')))
    fig.add_trace(go.Scatter(x=[log10alp] * 2, y=[0, xy_max], showlegend=False, hoverinfo='skip', mode='lines',
                             line=dict(color='gray', dash='solid')))
    # plot lines that indicate FDR caps
    fig.add_trace(go.Scatter(x=[0, xy_max], y=[log10FDR_limit] * 2, showlegend=False, hoverinfo='skip', mode='lines',
                             line=dict(color='gray', dash='dash')))
    fig.add_trace(go.Scatter(x=[log10FDR_limit] * 2, y=[0, xy_max], showlegend=False, hoverinfo='skip', mode='lines',
                             line=dict(color='gray', dash='dash')))
    # plot scatter points for the 4 groups of genes
    for name_i, color_i in zip(['both', 'neither', f'{method1} only', f'{method2} only'],
                               [color_both, color_neither, color_only1, color_only2]):
        ind_i = df_plot.Significant == name_i
        X = df_plot.loc[ind_i, f"{method1}: -log10(FDR) capped"].to_numpy()
        Y = df_plot.loc[ind_i, f"{method2}: -log10(FDR) capped"].to_numpy()
        labels = df_plot.loc[df_plot.Significant == name_i, "GENE"].tolist()
        if name_i != 'neither':
            hoverinfo = 'name+text'
            # redefining hoverinfo
            X_lab = df_plot.loc[ind_i, f"{method1}: -log10(FDR)"].tolist()
            Y_lab = df_plot.loc[ind_i, f"{method2}: -log10(FDR)"].tolist()
            labels = [(f"({X_lab[i]:.2f}, {Y_lab[i]:.2f})<br>{l}") for i, l in enumerate(labels)]
            # find groups of points that overlap
            pts_overlap, pt_groups = find_overlapping_points(X, Y)
            pts_nooverlap = np.setdiff1d(np.arange(len(X)), pts_overlap)
            # updating hoverinfo for groups of points that overlap
            for ptg in pt_groups:
                newlab = '<br>'.join(np.array(labels)[ptg].tolist())
                for i in ptg:
                    labels[i] = newlab
        else:
            hoverinfo = 'all'
        # plot scatter points
        fig.add_trace(go.Scatter(
            x=X.tolist(),
            y=Y.tolist(),
            text=labels,
            mode='markers',
            name=name_i,
            marker=dict(color=color_i, opacity=opacity),
            hoverinfo=hoverinfo,
            xhoverformat='.2f',
            yhoverformat='.2f'
        ))
        # add annotations
        if name_i != 'neither':
            # plot annotations for non-overlapping points
            for i in pts_nooverlap:
                fig.add_annotation(
                    x=X[i],
                    y=Y[i],
                    text=labels[i].split('<br>')[-1],
                    showarrow=True,
                    ax=10,
                    ay=7.5,
                    font=dict(color=color_i),
                    textangle=0,
                    xanchor="left",
                    yanchor="middle"
                )
            # plot annotations for overlapping points
            for ptg in pt_groups:
                ay_start = -ay_gap * (ptg.shape[0] - 1) / 2
                for j, i in enumerate(ptg):
                    fig.add_annotation(
                        x=X[i],
                        y=Y[i],
                        text=labels[i].split('<br>')[2 * j + 1],
                        showarrow=True,
                        ax=10,
                        ay=ay_start + j * ay_gap,
                        font=dict(color=color_i),
                        textangle=0,
                        xanchor="left",
                        yanchor="middle"
                    )
    # format axes and style
    fig.update_layout(
        width=fig_size[0],
        height=fig_size[1],
        title=dict(
            text='{} vs. {}'.format(method1, method2),
            x=0.5,
            xanchor='center'
        ),
        font=dict(size=fs),
        legend=dict(
            title="Gene significant for:",
            x=1.15,  # Move the legend further to the right
            xanchor='left'
        ),
        template='plotly_white',
        xaxis=dict(title=f'{method1}: -log10(FDR)', range=[0, xy_max], dtick=tick_step, scaleanchor='y'),
        # Link the scales of x and y
        yaxis=dict(title=f'{method2}: -log10(FDR)', range=[0, xy_max], dtick=tick_step)
    )
    # collect significant gene names for the three categories with such genes
    sig_genes_dict = {}
    sig_genes_dict[f'{method1} and {method2}'] = df_plot.loc[df_plot.Significant == 'both', 'GENE'].to_numpy()
    is_nan = df_plot.Significant.isna()
    ind_method1 = np.logical_and(is_nan, df_plot[f'{method1}: -log10(FDR)'] >= log10alp)
    sig_genes_dict[f'{method1} only'] = df_plot.loc[
                                            df_plot.Significant == f'{method1} only', 'GENE'].to_numpy().tolist() + \
                                        df_plot.loc[ind_method1, 'GENE'].to_numpy().tolist()
    ind_method2 = np.logical_and(is_nan, df_plot[f'{method2}: -log10(FDR)'] >= log10alp)
    sig_genes_dict[f'{method2} only'] = df_plot.loc[
                                            df_plot.Significant == f'{method2} only', 'GENE'].to_numpy().tolist() + \
                                        df_plot.loc[ind_method2, 'GENE'].to_numpy().tolist()
    return fig, sig_genes_dict


def plot_table_comparison(df, sig_genes_dict):
    """
    Generates and displays formatted Plotly tables to visualize significant genes
    based on different methods.

    Parameters
    ----------
    df : pandas.DataFrame
        Input DataFrame containing gene-related information. Expected columns include:
        - 'GENE': Gene names.
        - 'RANK': Ranking of genes.
        - 'SIZE_coding': Coding region size.
        - Significance metrics such as 'PVAL', 'FDR', and method-specific columns prefixed with 'SIG_'.

    sig_genes_dict : dict
        Dictionary mapping methods or conditions (e.g., "Method A only") to lists of
        significant genes identified under those methods.

    Workflow
    --------
    For each key in `sig_genes_dict`:
    1. Filter `df` to include only rows corresponding to genes in the associated list.
    2. Sort the filtered DataFrame by the 'RANK' column.
    3. Annotate genes untested by other methods with an asterisk (*) if applicable.
    4. Format numeric columns (e.g., P-values, FDR) for readability.
    5. Replace missing or invalid values (e.g., 'nan') with 'NA'.
    6. Create a Plotly `Table` figure:
       - Include bold headers.
       - Apply alternating row colors for better readability.
       - Align the first column ('GENE') to the left; center-align all other columns.
    7. Add a title to indicate the current method and include annotations for untested genes.

    Formatting
    ----------
    - Colors for headers, rows, and lines are controlled via global variables:
      `headerColor`, `rowOddColor`, `rowEvenColor`, and `lineColor`.
    - Numeric values are reformatted using an external function `reformat_numbers`.

    Returns
    -------
    list of plotly.graph_objs.Figure
        A list of Plotly `Figure` objects, one for each method in `sig_genes_dict`.

    Notes
    -----
    - The function assumes numeric columns requiring formatting contain substrings
      like 'PVAL' or 'FDR' in their names.
    - Asterisks (*) annotate genes not tested by other methods, with a footnote added
      to the table as an annotation.

    Example
    -------
    >>> sig_genes_dict = {
            "Method A only": ["GENE1", "GENE2"],
            "Method B only": ["GENE3", "GENE4"]
        }
    >>> tables = plot_table_comparison(df, sig_genes_dict)
    >>> for table in tables:
            table.show()

    This will display separate tables for "Method A only" and "Method B only", showing
    their respective significant genes.
    """
    methods = [k.split(' ')[0] for k in list(sig_genes_dict.keys()) if 'only' in k]
    fig_tables = []
    for key in sig_genes_dict:
        df_plot = df.loc[df.GENE.isin(sig_genes_dict[key])].copy().sort_values('RANK')
        df_plot = df_plot[['GENE', 'RANK', 'CHROM', 'SIZE_coding'] + [c for c in df_plot.columns if len([m for m in methods if m in c]) > 0] + ['CGC', 'PANCAN']]
        method_other = [m for m in methods if m not in key]
        if len(method_other) > 0:
            table_annot = df_plot['SIG_' + method_other[0]].isna().sum() > 0
            df_plot.loc[df_plot['SIG_' + method_other[0]].isna(), 'GENE'] += '*'
        else:
            table_annot = False
        # formatting data in columns
        df_plot = format_cols(df_plot)
        # generate table figure
        fig_table = go.Figure(data=[go.Table(
            header=dict(values=['<b>' + col + '</b>' for col in df_plot.columns],
                        line_color=lineColor,
                        fill_color=headerColor,
                        align=['left'] + ['center'] * (len(df_plot.columns) - 1),
                        font=dict(color='white', size=12)
                        ),
            cells=dict(values=[df_plot[col].tolist() for col in df_plot.columns],
                       line_color=lineColor,
                       fill_color=[[rowOddColor if i % 2 == 0 else rowEvenColor for i in range(df_plot.shape[0])]],
                       align=['left'] + ['center'] * (len(df_plot.columns) - 1),
                       font=dict(color=lineColor, size=11),
                       )
        )])
        # table height
        height = 100 + len(df_plot) * 20
        # add title
        fig_table.update_layout(
            title=dict(
                text='Significant with ' + key + ':',
                font=dict(size=18),
                x=0.5,
                y=(height-1)/height,
                xanchor='center',
                yanchor='top',
            ),
            margin=dict(r=5, l=5, t=30, b=30),
            height=height
        )
        if table_annot:
            fig_table.update_layout(
                annotations=[
                    dict(
                        text="{{GENE}}*: Gene not tested by {}".format(method_other[0]),
                        x=0,
                        y=-15/(height-60),
                        xref="paper",
                        yref="paper",
                        showarrow=False,
                        align="left",
                        valign="top",
                        font=dict(size=12)
                    )
                ]
            )
        fig_tables.append(fig_table)
    return fig_tables

def plot_table_summary(df):
    """
    Generate a series of Plotly tables summarizing significant genes across different methods.

    This function identifies subsets of methods based on columns with "SIG_" in their names,
    generates tables for genes significant in each subset, and displays the tables with
    enhanced formatting and optional annotations.

    Parameters
    ----------
    df : pandas.DataFrame
        The input DataFrame containing gene-related information. Expected columns include:
        - Columns with the prefix "SIG_" to indicate significance for various methods.
        - Columns with "PVAL" or "FDR" for numerical reformatting.
        - Columns "RANK" and "SIZE_coding" for integer formatting.
        - A "GENE" column to annotate genes.

    Returns
    -------
    fig_tables : dict of plotly.graph_objects.Figure
        A dictionary where keys are method subsets and values are Plotly table figures.
    """

    # find all method types
    methods = [c.split('_')[-1] for c in df.columns if 'SIG' in c]
    # generate all subsets using chain and combinations
    subsets = chain.from_iterable(combinations(methods, r) for r in range(len(methods) + 1))
    # sort subsets by length in descending order
    subsets = sorted(subsets, key=len, reverse=True)[:-1]

    fig_tables = {}
    for i, s in enumerate(subsets):
        title_short = ', '.join(s[:-2]) + (', ' if len(s) > 2 else '') + ' and '.join(s[-2:]) + (
            ' only' if i > 0 else '')
        title = 'Significant with ' + title_short + ':'

        is_sig = np.all(df[['SIG_' + m for m in s]] == True, axis=1)
        if i == 0:
            is_keep = is_sig
            other_methods = []
        else:
            other_methods = np.setdiff1d(subsets[0], s).tolist()
            is_notsig = ~np.any(df[['SIG_' + m for m in other_methods]] == True, axis=1)
            is_keep = np.logical_and(is_sig, is_notsig)
        dfi = df.loc[is_keep].copy()
        if len(other_methods) > 0:
            table_annot = dfi[['SIG_' + m for m in other_methods]].isna().to_numpy().sum() > 0
            dfi.loc[np.any(dfi[['SIG_' + m for m in other_methods]].isna(), axis=1), 'GENE'] += '*'
        else:
            table_annot = False
        # formatting data in columns
        dfi = format_cols(dfi)
        # generate table figure
        fig_table = go.Figure(data=[go.Table(
            header=dict(values=['<b>' + col + '</b>' for col in dfi.columns],
                        line_color=lineColor,
                        fill_color=headerColor,
                        align=['left'] + ['center'] * (len(dfi.columns) - 1),
                        font=dict(color='white', size=12)
                        ),
            cells=dict(values=[dfi[col].tolist() for col in dfi.columns],
                       line_color=lineColor,
                       fill_color=[[rowOddColor if i % 2 == 0 else rowEvenColor for i in range(dfi.shape[0])]],
                       align=['left'] + ['center'] * (len(dfi.columns) - 1),
                       font=dict(color=lineColor, size=11),
                       )
        )])
        # table height
        height = 100 + len(dfi) * 20
        # add title
        fig_table.update_layout(
            title=dict(
                text=title,
                font=dict(size=18),
                x=0.5,
                y=(height - 1) / height,
                xanchor='center',
                yanchor='top',
            ),
            margin=dict(r=5, l=5, t=30, b=30),
            height=height
        )
        if table_annot:
            fig_table.update_layout(
                annotations=[
                    dict(
                        text="{GENE}*: Gene not tested by at least one of the other method(s)",
                        x=0,
                        y=-15 / (height - 60),
                        xref="paper",
                        yref="paper",
                        showarrow=False,
                        align="left",
                        valign="top",
                        font=dict(size=12)
                    )
                ]
            )
        fig_tables[title_short] = fig_table
    return fig_tables

#
# MAIN FUNCTION THAT GENERATES REPORTS

def generate_reports(
        path_mutsig,
        path_dndscv,
        path_dig,
        path_cgc,
        path_pancan,
        dir_output,
        prefix_output=None,
        alp=0.1
):
    # Load and Format Output Files of Statistical Methods
    # collect and process results from different statistical methods
    df = preprocess_results(path_mutsig, path_dndscv, path_dig, path_cgc, path_pancan, alp)
    # save the processed results to a TSV file
    # df.to_csv(dir_output + '/' + ('' if (prefix_output is None) else prefix_output + '_') + 'merged_results.tsv', sep='\t', index=False)

    # Generate HTML That Compares Pairs of Statistical Methods

    # html_comparison = """
    #     <!DOCTYPE html>
    #     <html lang="en">
    #     <head>
    #         <meta charset="UTF-8">
    #         <meta name="viewport" content="width=device-width, initial-scale=1.0">
    #         <title>Compare Tools</title>
    #         <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    #         <style>
    #             h1 {{
    #                 margin-bottom: 20px;
    #             }}
    #             label, select {{
    #                 font-size: 18px;
    #             }}
    #             /* Browser Reset */
    #             * {{
    #                 margin: 0;
    #                 padding: 0;
    #                 box-sizing: border-box;
    #             }}
        
    #             /* Centering Figure Container */
    #             .figure-container {{
    #                 display: flex;                /* Flexbox layout */
    #                 justify-content: center;      /* Horizontal centering */
    #                 align-items: center;          /* Vertical centering */
    #                 width: 100%;                  /* Full width container */
    #                 margin: 10px 0;               /* Spacing around the figure */
    #             }}
        
    #             #fig-plot {{
    #                 display: inline-block;        /* Ensure it only takes the required width */
    #                 justify-content: center;      /* Horizontal centering */
    #                 align-items: center;          /* Vertical centering */
    #                 max-width: 50%;               /* Optional: Prevent it from being too wide */
    #             }}
        
    #             /* Full-Width Table Container */
    #             .table-container {{
    #                 width: 100%;                  /* Table spans full width */
    #                 margin: 0px 0;               /* Spacing around tables */
    #             }}
        
    #             .plot {{
    #                 width: 100%;                  /* Full width for tables */
    #             }}
    #         </style>
    #     </head>
    #     <body>
    #         <h1>Comparison of Driver Discovery Tool Results for Coding Regions</h1>
        
    #         <label for="methods">Tools to compare:</label>
    #         <select id="methods" onchange="updatePlot()">
    #             {methods_options}
    #         </select>
        
    #         <!-- Centered Figure -->
    #         <div class="figure-container">
    #             <div id="fig-plot" class="plot"></div>
    #         </div>
        
    #         <!-- Full-Width Tables -->
    #         <div class="table-container">
    #             <div id="table1-plot" class="plot"></div>
    #         </div>
    #         <div class="table-container">
    #             <div id="table2-plot" class="plot"></div>
    #         </div>
    #         <div class="table-container">
    #             <div id="table3-plot" class="plot"></div>
    #         </div>
        
    #         <script>
    #             var plotData = {plot_data};
        
    #             function updatePlot() {{
    #                 var methodsTypeKey = document.getElementById("methods").value;
        
    #                 var data = plotData[methodsTypeKey];
        
    #                 // Update Figure
    #                 var figData = data.fig;
    #                 Plotly.react('fig-plot', figData);
        
    #                 // Update Tables
    #                 var table1Data = data.table1;
    #                 Plotly.react('table1-plot', table1Data);
        
    #                 var table2Data = data.table2;
    #                 Plotly.react('table2-plot', table2Data);
        
    #                 var table3Data = data.table3;
    #                 Plotly.react('table3-plot', table3Data);
    #             }}
        
    #             // Initial plot
    #             updatePlot();
    #         </script>
    #     </body>
    #     </html>
    #     """
    # generate the dropdown options
    # methods_options = "\n".join([f'<option value="{key}">{key}</option>' for key in methods.keys()])
    # prepare plot data for all dropdown options
    plot_data = {}
    for methods_key, methods_val in methods.items():
        # generate FDR comparison plot
        fig, sig_genes_dict = plot_fdr_comparison(df, methods_val[0], methods_val[1], alp)
        # generate table plot
        fig_table = plot_table_comparison(df, sig_genes_dict)
        # store plotly figures
        plot_data[methods_key] = {
            'fig': fig.to_dict(),
        }
        for i, fig_table_i in enumerate(fig_table):
            plot_data[methods_key][f'table{i+1}'] = fig_table_i.to_dict()


    # convert plot data to JSON-like structure
    # plot_data_json = json.dumps(plot_data)
    # # combine everything into the final HTML
    # html_comparison = html_comparison.format(
    #     methods_options=methods_options,
    #     plot_data=plot_data_json
    # )
    # # save to an HTML file
    # comparison_path_rel = ('' if (prefix_output is None) else prefix_output + '_') + 'merged_report_comparison.html'
    # comparison_path = dir_output + '/' + comparison_path_rel
    # with open(comparison_path, 'w') as f:
    #     f.write(html_comparison)

    # Generate HTML That Summarizes Statistical Methods

    # html_summary = """
    #     <!DOCTYPE html>
    #     <html lang="en">
    #     <head>
    #         <meta charset="UTF-8">
    #         <meta name="viewport" content="width=device-width, initial-scale=1.0">
    #         <title>Summary</title>
    #         <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    #         <style>
    #             h1 {{
    #                 margin-bottom: 20px;
    #             }}
    #             h2 {{
    #                 margin-bottom: 10px;
    #             }}
    #             label, select {{
    #                 font-size: 18px;
    #             }}
    #             * {{
    #                 margin: 0;
    #                 padding: 0;
    #                 box-sizing: border-box;
    #             }}
    #             .table-container {{
    #                 display: flex;                /* Flexbox layout */
    #                 justify-content: center;      /* Horizontal centering */
    #                 align-items: center;          /* Vertical centering */
    #                 width: 100%;                  /* Full width container */
    #                 margin: 0px 0;               /* Spacing around the figure */
    #             }}
    #             .plot {{
    #                 width: 100%;                  /* Full width for tables */
    #             }}
    #         </style>
    #     </head>
    #     <body>
    #         <h1>Summary of Driver Discovery Tool Results for Coding Regions</h1>
            
    #         <h2>Genes Ranked According to Combined FDR</h2>
            
    #         <div class="table-container">
    #             <div class="plot">{table_static}</div>
    #         </div>
            
    #         <h2>Concordance Across Tools</h2>
            
    #         <label for="methods">List genes significant with:</label>
    #         <select id="methods" onchange="updatePlot()">
    #             {methods_options}
    #         </select>
    #         <div class="table-container">
    #             <div id="table-plot" class="plot"></div>
    #         </div>

    #         <script>
    #             var plotData = {plot_data};

    #             function updatePlot() {{
    #                 var methodsTypeKey = document.getElementById("methods").value;

    #                 var data = plotData[methodsTypeKey];

    #                 // Update Dynamic Table
    #                 var tableData = data.table;
    #                 Plotly.react('table-plot', tableData);
    #             }}

    #             // Initial plot
    #             updatePlot();
    #         </script>
    #     </body>
    #     </html>
    #         """

    # Genrate dynamic table of HTML:
    fig_tables = plot_table_summary(df)
    # generate the dropdown options
    methods_options = "\n".join([f'<option value="{key}">{key}</option>' for key in fig_tables.keys()])
    # prepare plot data for all dropdown options
    plot_data = {}
    for key in fig_tables:
        plot_data[key] = {
            'table': fig_tables[key].to_dict(),
        }
    # convert plot data to JSON-like structure
    plot_data_json = json.dumps(plot_data)

    # Generate static table of HTML:
    # find significant genes w.r.t combined FDR
    ind_sig = df['FDR_min'] <= alp
    n_rows = max(n_rows_min, int(np.sum(ind_sig) * (1 + n_rows_buffer)))
    df_table = df.iloc[:n_rows].copy()
    # generate table figure
    df_table = format_cols(df_table)
    # making the significant rows bold
    for i in range(df_table.shape[0]):
        if ind_sig[i]:
            df_table.loc[i, :] = '<b>' + df_table.loc[i, :].astype(str) + '</b>'

    fig_table_static = go.Figure(data=[go.Table(
        header=dict(values=['<b>' + col + '</b>' for col in df_table.columns],
                    line_color='darkslategray',
                    fill_color=headerColor,
                    align=['left'] + ['center'] * (len(df_table.columns) - 1),
                    font=dict(color='white', size=12)
                    ),
        cells=dict(values=[df_table[col].tolist() for col in df_table.columns],
                   line_color='darkslategray',
                   fill_color=[[rowOddColor if i % 2 == 0 else rowEvenColor for i in range(df_table.shape[0])]],
                   align=['left'] + ['center'] * (len(df_table.columns) - 1),
                   font=dict(color='darkslategray', size=11),
                   # format=['html'] * len(df_plot.columns)  # Enable HTML formatting
                   )
    )])

    # # save static figure as HTML div
    # table_static_html = fig_table_static.to_html(full_html=False, include_plotlyjs='cdn')

    # # combine everything into the final HTML
    # html_summary = html_summary.format(
    #     methods_options=methods_options,
    #     plot_data=plot_data_json,
    #     table_static=table_static_html
    # )
    # # save to an HTML file
    # summary_path_rel = ('' if (prefix_output is None) else prefix_output + '_') + 'merged_report_summary.html'
    # summary_path = dir_output + '/' + summary_path_rel
    # with open(summary_path, 'w') as f:
    #     f.write(html_summary)

    # html_main_report = f"""
    #     <!DOCTYPE html>
    #     <html lang="en">
    #     <head>
    #         <meta charset="UTF-8">
    #         <meta name="viewport" content="width=device-width">
    #         <title>Driver Discovery Suit Report</title>
    #         <style>
    #             .report-section {{
    #                 display: none;
    #             }}
    #             .active {{
    #                 display: block;
    #             }}
    #             .navbar {{
    #                 overflow: hidden;
    #                 background-color: #333;
    #             }}
    #             .navbar a {{
    #                 float: left;
    #                 display: block;
    #                 color: #f2f2f2;
    #                 text-align: center;
    #                 padding: 14px 16px;
    #                 text-decoration: none;
    #             }}
    #             .navbar a:hover {{
    #                 background-color: #ddd;
    #                 color: black;
    #             }}
    #             .navbar a.active-link {{
    #                 background-color: white;
    #                 color: black;
    #             }}
    #             iframe {{
    #                 width: 100%;
    #                 height: 1000px;
    #                 border: none;
    #             }}
    #         </style>
    #     </head>
    #     <body>
    #         <div class="navbar">
    #             <a href="#" onclick="showReport('summary', '{summary_path_rel}')">Summary</a>
    #             <a href="#" onclick="showReport('comparison', '{comparison_path_rel}')">Comparison</a>
    #         </div>

    #         <div id="summary" class="report-section active">
    #             <iframe id="report-frame" src="{summary_path_rel}"></iframe>
    #         </div>

    #         <script>
    #             function showReport(reportId, reportUrl) {{
    #                 var iframe = document.getElementById('report-frame');
    #                 iframe.src = reportUrl;

    #                 var links = document.querySelectorAll('.navbar a');
    #                 links.forEach(link => link.classList.remove('active-link'));

    #                 var activeLink = document.querySelector(`.navbar a[onclick*="${{reportId}}"]`);
    #                 activeLink.classList.add('active-link');
    #             }}
    #         </script>
    #     </body>
    #     </html>
    #     """
    # # save to an HTML file
    # with open(dir_output + '/' + ('' if (prefix_output is None) else prefix_output + '_') + f'merged_report_main.html',
    #           'w') as f:
    #     f.write(html_main_report)
