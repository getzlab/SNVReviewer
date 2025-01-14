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
mut_type = {
    'Indels + SNVs': 'MUT',
    'Indels': 'INDEL',
    'SNVs': 'SNV'
}
scatterpoint_type = {
    "Uniform P-mid": "unif",
    "P-mid": "recalc"
}
# display_bounds_type = {
#     'Yes': True,
#     'No': False
# }
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
text_special = 'Sample-wise case does not exist for Indels and Indels + Nonsynonymous SNVs!'

# derived parameters
col_err_sig = ','.join(col_sig.split(',')[:-1]) + ', {})'.format(opac_err)
col_err_nonsig = ','.join(col_nonsig.split(',')[:-1]) + ', {})'.format(opac_err)

# dictionaries for the two dropdowns
burden_type = {
    'Total': 'BURDEN',
    'Sample-wise': 'BURDEN_SAMPLE',
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

# not used for the coding report results
# mutation_type = {
#     'Indels + SNVs': 'MUT',
#     'Indels': 'INDEL',
#     'SNVs': 'SNV'
# }

# not used for the coding region plot generation
# scatterpoint_type = {
#     "Uniform P-mid": "unif",
#     "P-mid": "recalc"
# }
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


# def generate_dig_report(
#         dfi, 
#         dig_output_type,
#         # path_to_coding_results,
#         # path_to_promoter_results,
#         # path_to_3utr_results,
#         # path_to_5utr_results,
#         dir_output,
#         prefix_output=None,
#         alp=0.1
# ):
#     # Driver gene lists
#     cgc_list_path = "gs://getzlab-workflows-reference_files-oa/hg19/dig/cancer_gene_census_2024_06_20.tsv"
#     pancan_list_path = "gs://getzlab-workflows-reference_files-oa/hg19/dig/pancanatlas_genes.tsv"

#     # Driver gene lists
#     cgc_list = pd.read_csv(cgc_list_path, sep='\t').to_numpy().flatten()
#     pancan_list = pd.read_csv(pancan_list_path, sep='\t').to_numpy().flatten()
#     # DIG results
#     # dig_outputs = {
#     #     'coding': path_to_coding_results,
#     #     'promoter': path_to_promoter_results,
#     #     '5utr': path_to_5utr_results,
#     #     '3utr': path_to_3utr_results
#     # }
#     result_types = list(dig_outputs.keys())
#     cols_mt_set = set()

#     for result_type in result_types:

#         dfi = pd.read_csv(dig_outputs[result_type], sep='\t')

#         if result_type == 'coding':
#             dfi = dfi.set_index('GENE')
#             df_comb = dfi[['CHROM']].copy()
#             col_obs = ['OBS_NONSYN', 'N_SAMP_NONSYN']
#             col_pi = ['Pi_NONSYN', 'Pi_NONSYN']
#             col_size = 'GENE_LENGTH'

#         else:
#             dfi['GENE'] = dfi.ELT.str.split('::', expand=True)[2]
#             dfi = dfi.sort_values(['GENE', 'OBS_SNV'])
#             dfi = dfi.set_index('GENE')
#             dfi = dfi.loc[~dfi.index.duplicated()]
#             col_obs = ['OBS_SNV', 'OBS_SAMPLES']
#             col_pi = ['Pi_SUM', 'Pi_SUM']
#             col_size = 'ELT_SIZE'
#         cols_mt = ['SNV', 'SNV_SAMPLE']

#         if 'EXP_INDEL' in dfi.columns:
#             col_obs += ['OBS_INDEL']
#             col_pi += ['Pi_INDEL']
#             cols_mt += ['INDEL', 'MUT']

#         else:
#             for key in list(mut_type.keys()):
#                 if 'indel' in key.lower():
#                     del mut_type[key]

#         cols_mt_set.update(cols_mt)

#         genes_in_coding = dfi.index[dfi.index.isin(df_comb.index)]
#         dfi_comp = dfi.loc[genes_in_coding]
#         df_comb.loc[genes_in_coding, 'SIZE_' + result_type] = dfi_comp[col_size].to_numpy().copy()

#         for j, mt in enumerate(cols_mt[:3]):

#             if mt == 'INDEL':
#                 pfx_at = '_INDEL'
#             else:
#                 pfx_at = ''
#             df_comb.loc[genes_in_coding, 'PVAL_' + result_type + '_' + mt + '_recalc'] = nb_pvalue_greater_midp(
#                 dfi_comp[col_obs[j]],
#                 dfi_comp['ALPHA' + pfx_at],
#                 1 / (dfi_comp['THETA' + pfx_at] * dfi_comp[col_pi[j]] + 1)
#             )
#             df_comb.loc[genes_in_coding, 'PVAL_' + result_type + '_' + mt + '_unif'] = nb_pvalue_uniform_midp(
#                 dfi_comp[col_obs[j]],
#                 dfi_comp['ALPHA' + pfx_at],
#                 1 / (dfi_comp['THETA' + pfx_at] * dfi_comp[col_pi[j]] + 1)
#             )
#             df_comb.loc[genes_in_coding, 'PVAL_' + result_type + '_' + mt + '_lower'] = nb_pvalue_lower(
#                 dfi_comp[col_obs[j]],
#                 dfi_comp['ALPHA' + pfx_at],
#                 1 / (dfi_comp['THETA' + pfx_at] * dfi_comp[col_pi[j]] + 1)
#             )
#             df_comb.loc[genes_in_coding, 'PVAL_' + result_type + '_' + mt + '_upper'] = nb_pvalue_upper(
#                 dfi_comp[col_obs[j]],
#                 dfi_comp['ALPHA' + pfx_at],
#                 1 / (dfi_comp['THETA' + pfx_at] * dfi_comp[col_pi[j]] + 1)
#             )
#         if 'EXP_INDEL' in dfi.columns:
#             # combining p-values with Fisher's method
#             col_i = 'PVAL_' + result_type + '_' + 'MUT'
#             df_comb[col_i + '_recalc'] = np.nan
#             df_comb[col_i + '_unif'] = np.nan
#             df_comb[col_i + '_lower'] = np.nan
#             df_comb[col_i + '_upper'] = np.nan
#             for idx in df_comb.index:
#                 df_comb.at[idx, col_i + '_recalc'] = sp.stats.combine_pvalues(
#                     [df_comb.at[idx, 'PVAL_' + result_type + '_SNV_recalc'],
#                      df_comb.at[idx, 'PVAL_' + result_type + '_INDEL_recalc']], method='fisher')[1]
#                 df_comb.at[idx, col_i + '_unif'] = sp.stats.combine_pvalues(
#                     [df_comb.at[idx, 'PVAL_' + result_type + '_SNV_unif'],
#                      df_comb.at[idx, 'PVAL_' + result_type + '_INDEL_unif']], method='fisher')[1]
#                 df_comb.at[idx, col_i + '_lower'] = sp.stats.combine_pvalues(
#                     [df_comb.at[idx, 'PVAL_' + result_type + '_SNV_lower'],
#                      df_comb.at[idx, 'PVAL_' + result_type + '_INDEL_lower']], method='fisher')[1]
#                 df_comb.at[idx, col_i + '_upper'] = sp.stats.combine_pvalues(
#                     [df_comb.at[idx, 'PVAL_' + result_type + '_SNV_upper'],
#                      df_comb.at[idx, 'PVAL_' + result_type + '_INDEL_upper']], method='fisher')[1]

#     # combining p-values across region types
#     for mt in list(cols_mt_set):
#         for typ in ['recalc', 'unif', 'lower', 'upper']:
#             for idx in df_comb.index:
#                 df_comb.at[idx, 'PVAL_' + mt + '_' + typ] = sp.stats.combine_pvalues(
#                     [df_comb.at[idx, 'PVAL_' + rt + '_' + mt + '_' + typ] for rt in result_types if
#                      ~np.isnan(df_comb.at[idx, 'PVAL_' + rt + '_' + mt + '_' + typ])],
#                     method='fisher')[1]
#             isna = df_comb['PVAL_' + mt + '_' + typ].isna()
#             df_comb.loc[~isna, 'FDR_' + mt + '_' + typ] = fdrcorrection(df_comb.loc[~isna, 'PVAL_' + mt + '_' + typ])[1]

#     # adding CGC and PanCanAtlas information
#     df_comb['CGC'] = df_comb.index.isin(cgc_list).copy()
#     df_comb['PANCAN'] = df_comb.index.isin(pancan_list).copy()
#     df_comb = df_comb.reset_index().copy()

#     # CREATES THE COMBINED DIG REPORT RESULTS TXT FILE
#     # save combined p-values into a text file
#     df_comb.to_csv(('' if (prefix_output is None) else prefix_output + '.') + 'combined.dig.results.txt', sep='\t', index=False)

#     # generate the dropdown options
#     mut_options = "\n".join([f'<option value="{key}">{key}</option>' for key in mut_type.keys()])
#     burden_options = "\n".join([f'<option value="{key}">{key}</option>' for key in burden_type.keys()])
#     scatterpoint_options = "\n".join([f'<option value="{key}">{key}</option>' for key in scatterpoint_type.keys()])

#     # prepare plot data for all combinations of mut_type and burden_type dropdown options
#     plot_data = {}
#     for mut_key, mut_val in mut_type.items():
#         for bur_key, bur_val in burden_type.items():
#             for display_bounds_key, display_bounds_val in display_bounds_type.items():
#                 for scatterpoint_key, scatterpoint_val in scatterpoint_type.items():
#                     for display_labels_key, display_labels_val in display_labels_type.items():
#                         if not (mut_key in ['Indels', 'Indels + SNVs'] and bur_key == 'Sample-wise'):
#                             _, pvals, pval_bounds, labels, test_dom, ind_sig, table_fig = generate_plot_data(mut_val, bur_val, display_bounds_val, scatterpoint_val)

#                             # Q-Q Plot

#                             # Scatter plots
#                             x = -np.log10(np.arange(1, len(pvals) + 1) / (len(pvals) + 1))
#                             y = -np.log10(pvals)
#                             ind_capped = y > ymax
#                             ind_ncapped = np.logical_and(ind_sig, ~ind_capped)
#                             ind_nonsig = ~ind_sig
#                             ylim_upper = min(np.max(y), ymax) * (1 + hor_buffer)

#                             qq_fig = go.Figure()
#                             for i in range(len(result_types)):
#                                 ind = test_dom == result_types[i]
#                                 marker_i = markers[result_types[i]][1]
#                                 i_nonsig = np.logical_and(ind, ind_nonsig)
#                                 i_ncapped = np.logical_and(ind, ind_ncapped)
#                                 i_capped = np.logical_and(ind, ind_capped)

#                                 labels_capped = labels[i_capped].tolist()
#                                 x_capped = x[i_capped].tolist()
#                                 y_capped = y[i_capped].tolist()
#                                 x_ncapped = x[i_ncapped].tolist()
#                                 y_ncapped = y[i_ncapped].tolist()
#                                 x_nonsig = x[i_nonsig].tolist()
#                                 y_nonsig = y[i_nonsig].tolist()
#                                 if display_bounds_val:
#                                     y_upper = -np.log10(pval_bounds[:, 0]) - y
#                                     y_lower = y + np.log10(pval_bounds[:, 1])
#                                     y_upper_ncapped = y_upper[i_ncapped].tolist()
#                                     y_lower_ncapped = y_lower[i_ncapped].tolist()
#                                     y_upper_nonsig = y_upper[i_nonsig].tolist()
#                                     y_lower_nonsig = y_lower[i_nonsig].tolist()
#                                     dict_erry_sig = dict(
#                                         type='data',
#                                         symmetric=False,
#                                         array=y_upper_ncapped,
#                                         arrayminus=y_lower_ncapped,
#                                         thickness=thk_err,
#                                         width=wid_err,
#                                         color=col_err_sig
#                                     )
#                                     dict_erry_nonsig = dict(
#                                         type='data',
#                                         symmetric=False,
#                                         array=y_upper_nonsig,
#                                         arrayminus=y_lower_nonsig,
#                                         thickness=thk_err,
#                                         width=wid_err,
#                                         color=col_err_nonsig
#                                     )
#                                     labels_capped = get_labels_w_err(x_capped, y_capped, y_lower[i_capped].tolist(),
#                                                                      y_upper[i_capped].tolist(), labels_capped)
#                                     labels_ncapped = get_labels_w_err(x_ncapped, y_ncapped, y_lower_ncapped,
#                                                                      y_upper_ncapped, labels[i_ncapped].tolist())
#                                     labels_nonsig = get_labels_w_err(x_nonsig, y_nonsig, y_lower_nonsig, y_upper_nonsig,
#                                                                    labels[i_nonsig].tolist())
#                                     ylim_upper = min(np.max(y_upper + y), ymax) * (1 + hor_buffer)
#                                 else:
#                                     dict_erry_sig, dict_erry_nonsig = None, None
#                                     labels_capped = get_labels(x_capped, y_capped, labels_capped)
#                                     labels_ncapped = get_labels(x_ncapped, y_ncapped, labels[i_ncapped].tolist())
#                                     labels_nonsig = get_labels(x_nonsig, y_nonsig, labels[i_nonsig].tolist())
#                                     ylim_upper = min(np.max(y), ymax) * (1 + hor_buffer)

#                                 qq_fig = plot_qq(x_nonsig, y_nonsig, dict_erry_nonsig, labels_nonsig, col_nonsig, opac_nonsig,
#                                                  marker_i, 'Non-significant', qq_fig)
#                                 qq_fig = plot_qq(x_ncapped, y_ncapped, dict_erry_sig, labels_ncapped, col_sig, opac_sig,
#                                                  marker_i, 'Significant', qq_fig)
#                                 qq_fig = plot_qq(x_capped, [ymax] * sum(i_capped), None, labels_capped, col_sig,
#                                                  opac_sig, marker_i, 'Significant', qq_fig)

#                             # Dummy points for legend
#                             for test_i in markers.keys():
#                                 qq_fig.add_trace(
#                                     go.Scatter(
#                                         y=[None],
#                                         mode='markers',
#                                         marker=dict(
#                                             color='white',
#                                             symbol=markers[test_i][1],
#                                             size=msize * 1.25,
#                                             line=dict(color='black', width=2)
#                                         ),
#                                         name=markers[test_i][0]
#                                     )
#                                 )

#                             if display_labels_val:
#                                 x_capped = x[ind_capped].tolist()
#                                 x_ncapped = x[ind_ncapped].tolist()
#                                 y_capped = [ymax] * sum(ind_capped)
#                                 y_ncapped = y[ind_ncapped].tolist()
#                                 labels_capped = labels[ind_capped].tolist()
#                                 labels_ncapped = labels[ind_ncapped].tolist()
#                                 for (xi, yi, label) in zip(x_capped + x_ncapped, y_capped + y_ncapped,
#                                                            labels_capped + labels_ncapped):
#                                     qq_fig.add_annotation(
#                                         x=xi,
#                                         y=yi - ylim_upper * y_gap_annot,
#                                         text=label.split('<br>')[-1],
#                                         showarrow=False,
#                                         font=dict(color=col_sig),
#                                         textangle=-90,
#                                         xanchor="center",
#                                         yanchor="top"
#                                     )

#                             # Line plots
#                             qq_fig.add_trace(
#                                 go.Scatter(
#                                     x=[0, np.max(x)],
#                                     y=[0, np.max(x)],
#                                     mode='lines',
#                                     line=dict(dash=typ_thick, color=col_thick, width=thk_thick),
#                                     showlegend=False,
#                                     hoverinfo='skip'
#                                 )
#                             )
#                             if ymax < ylim_upper:
#                                 qq_fig.add_trace(
#                                     go.Scatter(
#                                         x=[0, np.max(x) * (1 + hor_buffer)],
#                                         y=[ymax] * 2,
#                                         mode='lines',
#                                         line=dict(dash=typ_thin, color=col_thin, width=thk_thin),
#                                         showlegend=False,
#                                         hoverinfo='skip'
#                                     )
#                                 )

#                             # Area of confidence intervals for the identity line
#                             xi = np.arange(1, len(pvals) + 1)
#                             clower = -np.log10(sp.stats.beta.ppf((1 - ci) / 2, xi, xi[::-1]))
#                             cupper = -np.log10(sp.stats.beta.ppf((1 + ci) / 2, xi, xi[::-1]))
#                             qq_fig.add_trace(go.Scatter(
#                                 x=np.concatenate([x, x[::-1]]).tolist(),  # Combine x values for fill
#                                 y=np.concatenate([clower, cupper[::-1]]).tolist(),  # Combine y values for fill
#                                 fill='toself',
#                                 fillcolor=col_beta,
#                                 line=dict(color='rgba(255,255,255,0)'),  # No line for the filled area
#                                 showlegend=False,
#                                 hoverinfo='skip'
#                             ))

#                             # Formatting the figure
#                             qq_fig.update_layout(
#                                 title='QQ-Plot of P-values:',
#                                 xaxis_title='Expected -Log10(P-value)',
#                                 yaxis_title='Observed -Log10(P-value)',
#                                 xaxis=dict(range=[0, np.max(x) * (1 + hor_buffer)]),
#                                 yaxis=dict(range=[0, ylim_upper]),
#                                 template='plotly_white',
#                                 legend=dict(
#                                     title=dict(
#                                         text="Dominant Region:",
#                                     ),
#                                     indentation=10
#                                 )
#                             )

#                             # Save figures as separate data
#                             plot_data[f"{mut_key}_{bur_key}_{display_bounds_key}_{scatterpoint_key}_{display_labels_key}"] = {
#                                 'qq': qq_fig.to_dict(),
#                                 'table': table_fig.to_dict(),
#                                 'text': bur_key + ' Mutation Burden of ' + mut_key,
#                                 'textcolor': 'black-text'
#                             }
#                         else:
#                             plot_data[f"{mut_key}_{bur_key}_{display_bounds_key}_{scatterpoint_key}_{display_labels_key}"] = {
#                                 'qq': None,
#                                 'table': None,
#                                 'text': text_special,
#                                 'textcolor': 'red-text'
#                             }

#     # convert plot data to JSON-like structure
#     plot_data_json = json.dumps(plot_data)

#     # combine everything into the final HTML
#     html_content = html_content.format(
#         mut_options=mut_options,
#         burden_options=burden_options,
#         scatterpoint_options=scatterpoint_options,
#         plot_data=plot_data_json
#     )

#     # save to an HTML file
#     with open(dir_output + '/' + ('' if (prefix_output is None) else prefix_output + '_') + f'dig_report_combined.html',
#               'w') as f:
#         f.write(html_content)



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

def generate_dig_report_plots(df):

    # prepare plot data for all combinations of mut_type and burden_type dropdown options
    plot_data = {}
    for mut_key, mut_val in mut_type.items():

        for bur_key, bur_val in burden_plot_type.items():

            for display_bounds_key, display_bounds_val in display_bounds_type.items():

                for scatterpoint_key, scatterpoint_val in scatterpoint_type.items():
                    for display_labels_key, display_labels_val in display_labels_type.items():
                        if not (mut_key in ['Indels', 'Indels + SNVs'] and bur_key == 'Sample-wise'):
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
                            plot_data[f"{mut_key}_{bur_key}_{display_bounds_key}_{scatterpoint_key}_{display_labels_key}"] = {
                                'qq': None,
                                'table': None,
                                'text': text_special,
                                'textcolor': 'red-text'
                            }
    
    return qq_fig, table_fig


# def parse_args():
#     """
#     Parse command-line arguments.
#     """
#     parser = argparse.ArgumentParser(description="Generate combined DIG report.")
#     parser.add_argument("path_to_coding_results", type=str, help="Path to the DIG results for coding regions.")
#     parser.add_argument("path_to_promoter_results", type=str, help="Path to the DIG results for promoter regions.")
#     parser.add_argument("path_to_3utr_results", type=str, help="Path to the DIG results for 3' UTR regions.")
#     parser.add_argument("path_to_5utr_results", type=str, help="Path to the DIG results for 5' UTR regions.")
#     parser.add_argument("dir_output", type=str, help="Output directory.")
#     parser.add_argument("cgc_list", type=str, help="Path to the list of CGC genes.")
#     parser.add_argument("pancan_list", type=str, help="Path to the list of PanCanAtlas genes.")
#     parser.add_argument("--prefix_output", type=str, default=None, help="Prefix for the output file.")
#     parser.add_argument("--alp", type=float, default=0.1, help="Significance level (default: 0.1).")
#     return parser.parse_args()
