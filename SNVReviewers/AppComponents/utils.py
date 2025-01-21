import argparse
import pandas as pd
import plotly.graph_objects as go
import numpy as np
import scipy as sp
from statsmodels.stats.multitest import fdrcorrection
import json

import pandas as pd
import numpy as np

# import plotly.graph_objects as go
import numpy as np
import scipy as sp
# from statsmodels.stats.multitest import fdrcorrection

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

# # properties of significant points
# col_sig = 'rgba(255, 0, 0, 1)'
# opac_sig = 0.7
# # properties of non-significant points
# col_nonsig = 'rgba(0, 0, 0, 1)'
# opac_nonsig = 0.5
# # properties of thinner lines
# thk_thin = 0.75
# col_thin = 'gray'
# typ_thin = 'dash'
# # properties of thicker lines
# thk_thick = 2.5
# col_thick = 'gray'
# typ_thick = 'dash'
# # properties of error bars
# wid_err = 2
# thk_err = 0.3
# opac_err = 0.3
# # properties of the bar plot
# col_bar = 'gray'
# opac_bar = 0.8

# color of beta confidence area
col_beta = 'rgba(128, 128, 128, 0.2)'

# # text for the special case when Sample-wise case does not exist
# text_special = 'Sample-wise case does not exist for Indels and Indels + SNVs!'

# # derived parameters
# col_err_sig = ','.join(col_sig.split(',')[:-1]) + ', {})'.format(opac_err)
# col_err_nonsig = ','.join(col_nonsig.split(',')[:-1]) + ', {})'.format(opac_err)

# markers used to indicate the dominant interval set in the QQ plot
markers = {
    'coding': ["Coding", "circle"],
    'promoter': ["Promoter", "square"],
    '5utr': ["5\' UTR", "star"],
    '3utr': ["3\' UTR", "triangle-down"]
}
msize = 7.5

# dropdown options
burden_plot_type = {
    'Total': '',
    'Sample-wise': 'SAMPLE',
}
result_types = ['coding', 'promoter', '5utr', '3utr']

scatterpoint_type = {
    "uniform_p_mid": "unif",
    "p_mid": "recalc"
}

display_labels_type = {
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

# properties of significant points
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



# dictionaries for the two dropdowns
burden_type = {
    'total': 'BURDEN',
    'sample_wise': 'BURDEN_SAMPLE',
}

mutation_type = {
    'Indels + Nonsynonymous SNVs': 'MUT',
    'Indels': 'INDEL',
    'Nonsynonymous SNVs': 'NONSYN',
    'Missense SNVs': 'MIS',
    'Nonsense SNVs': 'NONS',
    'Truncating SNVs': 'TRUNC',
    'Splice site SNVs': 'SPL',
    'Synonymous SNVs': 'SYN',
}

# used for generating the Combined QQ plot
mut_type = {
    'indels_snvs': 'MUT',
    'indels': 'INDEL',
    'snvs': 'SNV'
}

display_bounds_type = {
    'Yes': True,
    'No': False
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

def reformat_numbers(x, format='{:.3E}'):
    """
    Reformat numbers in an array to a specific format
    """
    return [format.format(n) for n in x]


# DIG REPORT FUNCTIONS
def generate_dig_report_dataframe(path_to_dig_results, alp=0.1):
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
    muts_ts = list(mutation_type.values())

    if 'EXP_INDEL' in df.columns:

        # Adding new columns for Non-synonymous SNVs + Indels
        df['OBS_MUT'] = df['OBS_NONSYN'] + df['OBS_INDEL']
        df['EXP_MUT'] = df['EXP_NONSYN'] + df['EXP_INDEL']
    else:
        for key in list(mutation_type.keys()):
            if 'indel' in key.lower():
                del mutation_type[key]

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


def reformat_numbers(df, cols, form='{:.2E}'):
    """
    Reformat numbers in an array to a specific format
    """
    return df[cols].applymap(lambda x: form.format(x) if not pd.isna(x) else 'NA')

# Make all the markers the same, i.e. don't have the dominant region, make sure to redownload the code from the google drive
def generate_plot_data(df, mut, bur, display_bounds, scatterpoint):
    """
    Given a mutation type and a burden type, generate the data for the volcano plot, Q-Q plot, and table plot
    :param mut: str, mutation type
    :param bur: str, burden type
    :param display_bounds: bool, whether to display the bounds of the p-values
    :param scatterpoint: str, type of p-values to use
    """
    df_comb = df.copy()
    alp = 0.1 # significance level
    # result_types = ['coding', 'promoter', '5utr', '3utr']

    col_chosen = mut + ('_SAMPLE' if bur=='SAMPLE' else '') + '_' + scatterpoint
    col_pvals = ['PVAL_' + rt + '_' + col_chosen for rt in result_types]
    col_sizes = ['SIZE_' + rt  for rt in result_types]
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
    df_plot.rename(columns={'PVAL_' + rt + '_' + col_chosen: 'PVAL_' + rt for rt in result_types}, inplace=True)

    df_plot[col_sizes] = reformat_numbers(df_plot, col_sizes, form='{:.0f}')
    cols_floats = ['PVAL', 'FDR'] + ['PVAL_' + rt for rt in result_types]
    df_plot[cols_floats] = reformat_numbers(df_plot, cols_floats)

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
def generate_dig_report_plots(df, mut_key, bur_key, display_bounds_key, display_labels_key, scatterpoint_key):

    # prepare plot data for all combinations of mut_type and burden_type dropdown options
    plot_data = {}
    
    mut_val = mut_type[mut_key]
    bur_val = burden_type[bur_key]
    display_bounds_val = display_bounds_type[display_bounds_key]
    scatterpoint_val = scatterpoint_type[scatterpoint_key]
    display_labels_val = display_labels_type[display_labels_key]
    text_special = ""

    if not (mut_key in ['indels', 'indels_snvs'] and bur_key == 'sample_wise'):
        _, pvals, pval_bounds, labels, test_dom, ind_sig, table_fig = generate_plot_data(df, mut_val, bur_val, display_bounds_val, scatterpoint_val)

        # Q-Q Plot
        # Scatter plots
        x = -np.log10(np.arange(1, len(pvals) + 1) / (len(pvals) + 1))
        y = -np.log10(pvals)
        ind_capped = y > ymax
        ind_ncapped = np.logical_and(ind_sig, ~ind_capped)
        ind_nonsig = ~ind_sig
        ylim_upper = min(np.max(y), ymax) * (1 + hor_buffer)

        qq_fig = go.Figure()

        for i in range(len(result_types)):

            ind = test_dom == result_types[i]
            marker_i = markers[result_types[i]][1]
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

        # CHECK WITH DAVID IF HE WANTS THE MARKERS ON THE FINAL VERSION OF THE QQ PLOT
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


# def reformat_numbers(x, format='{:.3E}'):
#     """
#     Reformat numbers in an array to a specific format
#     """
#     return [format.format(n) for n in x]


def generate_dig_report_coding(path_to_dig_results, dir_output, cgc_list_path, pancan_list_path, prefix_output=None, alp=0.1):
    # Driver gene lists
    cgc_list = pd.read_csv(cgc_list_path, sep='\t').to_numpy().flatten()
    pancan_list = pd.read_csv(pancan_list_path, sep='\t').to_numpy().flatten()
    # Output from DIGDriver
    df = pd.read_csv(path_to_dig_results, sep='\t')
    # df = df.iloc[:20]
    # Adding indicator of genes being part of the CGC or PanCan list
    df['CGC'] = df['GENE'].isin(cgc_list)
    df['PANCAN'] = df['GENE'].isin(pancan_list)
    muts_ts = list(mut_type.values())
    if 'EXP_INDEL' in df.columns:
        # Adding new columns for Non-synonymous SNVs + Indels
        df['OBS_MUT'] = df['OBS_NONSYN'] + df['OBS_INDEL']
        df['EXP_MUT'] = df['EXP_NONSYN'] + df['EXP_INDEL']
    else:
        for key in list(mut_type.keys()):
            if 'indel' in key.lower():
                del mut_type[key]
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

    def generate_plot_data_coding(mut, bur, display_bounds, scatterpoint):
        """
        Given a mutation type and a burden type, generate the data for the volcano plot, Q-Q plot, and table plot
        :param mut: str, mutation type
        :param bur: str, burden type
        :param display_bounds: bool, whether to display the bounds of the p-values
        :param scatterpoint: str, type of p-values to use
        """

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
            df_plot[col] = reformat_numbers(df_plot[col].to_numpy())
        for col in ['MU', 'SIGMA']:
            df_plot[col] = reformat_numbers(df_plot[col].to_numpy(), format='{:.2f}')
        for col in ['dNdS_EXP', 'EXP']:
            df_plot[col] = reformat_numbers(df_plot[col].to_numpy(), format='{:.3f}')

        df_plot['OBS'] = df_plot['OBS'].astype(int)
        is_inf = np.logical_or(np.isinf(df_plot.dNdS_OBS.to_numpy()), np.isnan(df_plot.dNdS_OBS.to_numpy()))
        dnds_obs = reformat_numbers(df_plot.loc[~is_inf, 'dNdS_OBS'].to_numpy().copy(), format='{:.3f}')
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
                       # format=['html'] * len(df_plot.columns)  # Enable HTML formatting
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

    # generate the data and table for the default values
    df_kept, pvals, pval_bounds, logfc, logq, logq_bounds, labels, ind_kept, table_fig = generate_plot_data('MIS',
                                                                                                            'BURDEN_SAMPLE',
                                                                                                            True,
                                                                                                            'recalc')

    # Update the HTML template to include the additional text below the table
    html_content = """
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>DIG for Coding Regions</title>
        <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
        <style>
            .container {{
                display: flex;
            }}
            .plot {{
                margin: 10px;
            }}
            #dnds-plot {{
                width: 33%;
            }}
            #volcano-plot {{
                width: 50%;
            }}
            #qq-plot {{
                width: 50%;
            }}
            #table-plot {{
                width: 100%;
            }}
            .figures-container {{
                display: flex;
                justify-content: space-between;
            }}
            .figure-column {{
                width: 33%;
            }}
            .switch {{
            position: relative;
            display: inline-block;
            width: 30px;
            height: 17px;
            }}
            .switch input {{
                opacity: 0;
                width: 0;
                height: 0;
            }}
            .slider {{
                position: absolute;
                cursor: pointer;
                top: 0;
                left: 0;
                right: 0;
                bottom: 0;
                background-color: #ccc;
                transition: .4s;
            }}
            .slider:before {{
                position: absolute;
                content: "";
                height: 13px;
                width: 13px;
                left: 2px;
                bottom: 2px;
                background-color: white;
                transition: .4s;
            }}
            input:checked + .slider {{
                background-color: #000000;
            }}
            input:checked + .slider:before {{
                transform: translateX(13px);
            }}
            .slider.round {{
                border-radius: 17px;
            }}
            .slider.round:before {{
                border-radius: 50%;
            }}
            .red-text {{
                color: red;
            }}
            .black-text {{
                color: black;
            }}
        </style>
    </head>
    <body>
        <h1>Mutation Burden Test for Coding Regions</h1>

        <label for="mut-type">Select Mutation Type:</label>
        <select id="mut-type" onchange="updatePlot()">
            {mut_options}
        </select>

        <label for="burden-type">Select Burden Type:</label>
        <select id="burden-type" onchange="updatePlot()">
            {burden_options}
        </select>

        <label for="scatterpoint-type">P-value Type:</label>
        <select id="scatterpoint-type" onchange="updatePlot()">
            {scatterpoint_options}
        </select>

        <label for="display-bounds-type">Display bounds:</label>
        <label class="switch">
            <input type="checkbox" id="display-bounds-type" onchange="updatePlot()">
            <span class="slider round"></span>
        </label>

        <h2 id="plot-title" class="text-color"></h2>

        <div class="container">
            <div id="volcano-plot" class="plot"></div>
            <div id="qq-plot" class="plot"></div>
        </div>
        <div id="table-plot" class="plot"></div>
        <div class="figures-container">
            <div class="figure-column">{fig_mu_html}</div>
            <div class="figure-column">{fig_sigma_html}</div>
            <div id="dnds-plot" class="plot"></div>
        </div>

        <script>
            var plotData = {plot_data};

            function updatePlot() {{
                var mutTypeKey = document.getElementById("mut-type").value;
                var burdenTypeKey = document.getElementById("burden-type").value;
                var scatterpointTypeKey = document.getElementById("scatterpoint-type").value;
                var displayBoundsKey = document.getElementById("display-bounds-type").checked ? 'Yes' : 'No';

                var data = plotData[mutTypeKey + '_' + burdenTypeKey + '_' + displayBoundsKey + '_' + scatterpointTypeKey];

                // Update Volcano Plot
                var volcanoData = data.volcano;
                Plotly.react('volcano-plot', volcanoData);

                // Update Q-Q Plot
                var qqData = data.qq;
                Plotly.react('qq-plot', qqData);

                // Update dNdS Plot
                var dndsData = data.dnds;
                Plotly.react('dnds-plot', dndsData);

                // Update Table Plot
                var tableData = data.table;
                Plotly.react('table-plot', tableData);

                // Update Header
                var headerText = data.text;
                document.getElementById("plot-title").textContent = headerText;
                var headerColor = data.textcolor;
                document.getElementById("plot-title").className = headerColor;
            }}

            // Initial plot
            updatePlot();
        </script>
    </body>
    </html>
    """

    # generate the dropdown options
    mut_options = "\n".join([f'<option value="{key}">{key}</option>' for key in mut_type.keys()])
    burden_options = "\n".join([f'<option value="{key}">{key}</option>' for key in burden_type.keys()])
    display_bounds_options = "\n".join([f'<option value="{key}">{key}</option>' for key in display_bounds_type.keys()])
    scatterpoint_options = "\n".join([f'<option value="{key}">{key}</option>' for key in scatterpoint_type.keys()])

    # prepare plot data for all combinations of mut_type and burden_type dropdown options
    plot_data = {}
    for mut_key, mut_val in mut_type.items():
        for bur_key, bur_val in burden_type.items():
            for display_bounds_key, display_bounds_val in display_bounds_type.items():
                for scatterpoint_key, scatterpoint_val in scatterpoint_type.items():
                    if not (mut_key in ['Indels', 'Indels + Nonsynonymous SNVs'] and bur_key == 'Sample-wise'):
                        df_kept, pvals, pval_bounds, logfc, logq, logq_bounds, labels, ind_kept, table_fig = generate_plot_data(
                            mut_val, bur_val, display_bounds_val, scatterpoint_val)

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
                        plot_data[f"{mut_key}_{bur_key}_{display_bounds_key}_{scatterpoint_key}"] = {
                            'volcano': None,
                            'qq': None,
                            'dnds': None,
                            'table': None,
                            'text': text_special,
                            'textcolor': 'red-text'
                        }

    # convert plot data to JSON-like structure
    plot_data_json = json.dumps(plot_data)

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

    # save static figures as HTML divs
    fig_mu_html = fig_mu.to_html(full_html=False, include_plotlyjs='cdn')
    fig_sigma_html = fig_sigma.to_html(full_html=False, include_plotlyjs='cdn')

    # combine everything into the final HTML
    html_content = html_content.format(
        mut_options=mut_options,
        burden_options=burden_options,
        display_bounds_options=display_bounds_options,
        scatterpoint_options=scatterpoint_options,
        plot_data=plot_data_json,
        fig_mu_html=fig_mu_html,
        fig_sigma_html=fig_sigma_html
    )

    # save to an HTML file
    with open(dir_output + '/' + ('' if (prefix_output is None) else prefix_output + '_') +
              'dig_report_coding_regions.html', 'w') as f:
        f.write(html_content)



