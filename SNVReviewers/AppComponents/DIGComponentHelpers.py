import plotly.graph_objects as go
import numpy as np
import pandas as pd
import scipy as sp

from SNVReviewers.AppComponents.utils import combined_burden_plot_type, combined_mutation_type, scatterpoint_type, display_bounds_type, combined_result_types, plot_qq
from SNVReviewers.AppComponents.utils import generate_coding_region_report, coding_region_mutation_type, coding_region_burden_type, markers, get_labels_w_err, get_labels
from SNVReviewers.AppComponents.utils import generate_dig_non_coding_plots, non_coding_region_burden_type, reformat_numbers, combined_reformat_numbers, msize, col_beta
from SNVReviewers.AppComponents.utils import nb_pvalue_lower, nb_pvalue_upper, nb_pvalue_uniform_midp, nb_pvalue_greater_midp, nb_pvalue_lower, nb_pvalue_upper

DIG_CODING_REGION_REPORT_COLUMN_NAMES = ["RANK", "GENE", 'CHROM', 'LENGTH', "PVAL", "FDR", "OBS", 
                                         "EXP", 'MU', 'SIGMA', 'dNdS_OBS', 'dNdS_EXP', 'FLAG', 'CGC', 'PANCAN']

DIG_REPORT_COMBINED_COLUMN_NAMES = ["RANK", "GENE", "FDR", "PVAL", "PVAL_coding", "PVAL_promoter", "PVAL_5utr", 
                                    "SIZE_coding", "SIZE_promoter", "SIZE_5utr", "SIZE_3utr", "CGC", "PANCAN"]

DIG_NON_CODING_REPORT_COLUMN_NAMES = ["RANK", "GENE", 'SIZE', "PVAL", "FDR", "OBS", 
                                      'EXP', 'MU', 'SIGMA', 'FLAG', 'CGC', 'PANCAN']

CODING_REGION_MUT_DROPDOWN = [
                        {'label': 'Indels + Nonsynonymous SNVs', 'value': 'indels_nonsynonymous_snvs'},
                        {'label': 'Indels', 'value': 'indels'},
                        {'label': 'Nonsynonymous + SNVs', 'value': 'nonsynonymous_snvs'},
                        {'label': 'Missense SNVs', 'value': 'missense_snvs'},
                        {'label': 'Nonsense SNVs', 'value': 'nonsense_snvs'},
                        {'label': 'Truncating SNVs', 'value': 'truncating_snvs'},
                        {'label': 'Splice Site SNVs', 'value': 'splice_site_snvs'},
                        {'label': 'Synonymous + SNVs', 'value': 'synonymous_snvs'}
                    ]
CODING_REGION_BUR_DROPDOWN = [
                        {'label': 'Total', 'value': 'total'},
                        {'label': 'Sample-wise', 'value': 'sample_wise'}
                    ]
CODING_REGION_SCATTER_DROPDOWN = [
                        {'label': 'Uniform P-mid', 'value': 'uniform_p_mid'},
                        {'label': 'P-mid', 'value': 'p_mid'},
                    ]

COMBINED_MUT_DROPDOWN = [
            {'label': 'Indels + SNVs', 'value': 'indels_snvs'},
            {'label': 'Indels', 'value': 'indels'},
            {'label': 'SNVs', 'value': 'snvs'}
        ]
COMBINED_BUR_DROPDOWN = [
            {'label': 'Total', 'value': 'total'},
            {'label': 'Sample-wise', 'value': 'sample_wise'}
        ]
COMBINED_SCATTER_DROPDOWN = [
            {'label': 'Uniform P-mid', 'value': 'uniform_p_mid'},
            {'label': 'P-mid', 'value': 'p_mid'},
        ]

SHOW_QQ_BIG_PLOT_STYLE = {'width':'1200px'}
SHOW_QQ_SMALL_PLOT_STYLE = {'width':'600px'}
HIDE_PLOT_STYLE = {'display':'none'}
SHOW_VOLCANO_PLOT_STYLE = {'display':'inline-block', 'width':'600px'}
SHOW_THREE_PLOTS_STYLE = {'display':'block', 'width':'350px'}
SHOW_TWO_PLOTS_STYLE = {'display':'block', 'width':'600'}

DIG_REPORT_VALUES = ["Combined", "Coding region", "Promoter region", "5-prime UTRs", "3-prime UTRs"]
DIG_DATAFRAME_IDX = 0
DND_DATAFRAME_IDX = 1
MUTSIG_DATAFRAME_IDX = 2
SNV_DATA_COLUMN_NAME = "snv_data"
DIG_LABEL_IDX = 1

# text for the special case when Sample-wise case does not exist
SPECIAL_TEXT = 'Sample-wise case does not exist for Indels and Indels + Nonsynonymous SNVs!'

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

# properties of near-significant points
col_nearsig = 'rgba(30, 144, 255, 1)'

# dropdown options
combined_burden_plot_type = {
    'total': '',
    'sample_wise': 'SAMPLE',
}

# derived parameters
col_err_sig = ','.join(col_sig.split(',')[:-1]) + ', {})'.format(opac_err)
col_err_nearsig = ','.join(col_nearsig.split(',')[:-1]) + ', {})'.format(opac_err)
col_err_nonsig = ','.join(col_nonsig.split(',')[:-1]) + ', {})'.format(opac_err)

# Make all the markers the same, i.e. don't have the dominant region, make sure to redownload the code from the google drive
def generate_combined_plot_data(df, mut, bur, display_bounds, scatterpoint):
    """
    Given a mutation type and a burden type, generate the data for the volcano plot, Q-Q plot, and table plot

    Parameters
    ==========
        df: pd.DataFrame
            DIG report data used to generate plots for AnnoMate dashboard
        mut: str, 
            mutation type
        bur: str, 
            burden type
        display_bounds: bool, 
            whether to display the bounds of the p-values
        scatterpoint: str, 
            type of p-values to use

    Return
    ======
        pd.DataFrame
            dataframe with the data to be displayed in the plotly table for the combined report type
        np.array 
            1D array of all the p values
        np.array 
            2D array of the upper and lower bounds for the p values
        np.array
            1D array of the gene names 
        np.array
            1D array with the dominant result type
        np.array
            1D array of the significant genes
        go.Figure
            plotly figure of the combined report type table
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

    # Ensures that rank, gene, cgc, and pancan columns are strings
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
                    )
    )])

    return df_kept, pvals, pval_bounds, labels, test_dom, ind_sig, table_fig

# for combined dig report plots
def generate_combined_dig_report_plots(
        df, 
        mut_key, 
        bur_key, 
        display_bounds_key, 
        scatterpoint_key
    ):
    """
    Generate the plots for the combined dig report type

    Parameters
    ==========
        df: pd.DataFrame,
            dataframe with the combined dig report data 
        mut_key: str,
            key that get used to access a specific value in the mutation dictionary 
        bur_key: str,
            key that get used to access a specific value in the burden dictionary  
        display_bounds_key: str,
            key that get used to access a specific value in the display bounds dictionary 
        scatterpoint_key: str,
            key that get used to access a specific value in the scatterpoint dictionary

    Return
    ======
        go.Figure
            plotly figure of the qq plot
        go.Figure
            plotly figure of the combined dig report table 
        str
            string with warning message if the user tries to view a bad mutation and burden type combination
    """
    mut_val = combined_mutation_type[mut_key]
    bur_val = combined_burden_plot_type[bur_key]
    display_bounds_val = display_bounds_type[display_bounds_key]
    scatterpoint_val = scatterpoint_type[scatterpoint_key]
    text_special = ""

    if not (mut_key in ['indels', 'indels_snvs'] and bur_key == 'sample_wise'):
        _, pvals, pval_bounds, labels, test_dom, ind_sig, table_fig = generate_combined_plot_data(df, mut_val, bur_val, display_bounds_val, scatterpoint_val)

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

    else:
        qq_fig = go.Figure()
        table_fig = go.Figure()
        text_special = SPECIAL_TEXT
    
    return qq_fig, table_fig, text_special

def gen_combined_app_component(
    dig_df,
    mutation_type,
    burden_type,
    p_val_type,
    display_bounds_value,
):
    """
    Generates all the content for the dig combined report type selection

    Parameters
    ==========
        dig_df: pd.DataFrame,
            dataframe with the dig report data
        mutation_type: str,
            string value of the mutation type that the user is viewing
        burden_type: str,
            string value of the burden type that the user is viewing
        p_val_type: str,
            string value of the p value type that the user is viewing
        display_bounds_value: str,
            string value of whether user wants to view the lower and upper bounds of the p values

    Return
    ======
        List[
            dig_df.to_dict('records'),
            volcano_fig,
            qq_fig,
            fig_mu, 
            fig_sigma,
            dnds_fig,
            dig_data_columns,
            mutation_type,
            burden_type,
            p_val_type,

            # updates the dropdown options
            COMBINED_MUT_DROPDOWN,
            COMBINED_BUR_DROPDOWN,
            COMBINED_SCATTER_DROPDOWN, 

            # changing the plot size and the visibility of the plots
            SHOW_QQ_BIG_PLOT_STYLE, # displays the qq plot
            HIDE_PLOT_STYLE,        # hides the volcano plot
            HIDE_PLOT_STYLE,        # hides the volcano plot
            HIDE_PLOT_STYLE,        # hides the volcano plot
            HIDE_PLOT_STYLE,        # hides the volcano plot

            # Warning text message (only generated if you get a specific mutation and burden type combination on dropdown menu)
            text_special,
            debugging
        ]   
    """
    debugging= ""
    all_page_content = []
    dig_df = dig_df.copy()
    dig_df = dig_df.sort_values(by='PVAL'+ "_"+ combined_mutation_type[mutation_type] + "_" + scatterpoint_type[p_val_type])
    dig_df['RANK'] = np.array([i+1 for i in range(len(dig_df))])

    display_bounds = display_bounds_value # whether to display the bounds
    
    qq_fig = go.Figure()
    fig_mu = go.Figure()
    fig_sigma = go.Figure()
    dnds_fig = go.Figure()
    volcano_fig = go.Figure()

    # checks if the user wants to display the bounds on the dig report plot
    if display_bounds_value:
        display_bounds = 'Yes'

    # defaults to not displaying lower/upper bounds on the dig report plot
    else:
        display_bounds = 'No'

    qq_fig, table_fig, text_special = generate_combined_dig_report_plots(dig_df, mutation_type, burden_type, display_bounds, p_val_type)
    dig_data_columns = []
    
    for clm_nm in DIG_REPORT_COMBINED_COLUMN_NAMES:
        dig_data_clm_dict = {}
        dig_data_clm_dict["name"] = clm_nm
        dig_data_clm_dict["id"] = clm_nm

        if clm_nm == 'FDR' or "PVAL" in clm_nm:

            # gets the column name corresponding to the mutation type and the scatterpoint type
            dig_data_clm_dict["id"] = clm_nm + "_"+ combined_mutation_type[mutation_type] + "_" + scatterpoint_type[p_val_type]

        dig_data_columns.append(dig_data_clm_dict)
    
    for column_dict in dig_data_columns:
        # gets the dig data column name
        column = column_dict["id"]

        # reformats the columns with float values in them
        if pd.api.types.is_float_dtype(dig_df[column]) and column != 'RANK':
            dig_df[column] = reformat_numbers(dig_df[column], format='{:.3E}')

    # ONLY GETTING THE FIRST 100 ROWS OF DATA TO DISPLAY IN THE TABLE
    # REMOVE THE DEBUGGING LATER!!!   
    dig_df = dig_df[:100]

    all_page_content = [
        dig_df.to_dict('records'),
        volcano_fig,
        qq_fig,
        fig_mu, 
        fig_sigma,
        dnds_fig,
        dig_data_columns,
        mutation_type,
        burden_type,
        p_val_type,

        # updates the dropdown options
        COMBINED_MUT_DROPDOWN,
        COMBINED_BUR_DROPDOWN,
        COMBINED_SCATTER_DROPDOWN, 

        # changing the plot size and the visibility of the plots
        SHOW_QQ_BIG_PLOT_STYLE, # displays the qq plot
        HIDE_PLOT_STYLE,        # hides the volcano plot
        HIDE_PLOT_STYLE,        # hides the volcano plot
        HIDE_PLOT_STYLE,        # hides the volcano plot
        HIDE_PLOT_STYLE,        # hides the volcano plot

        # Warning text message (only generated if you get a specific mutation and burden type combination on dropdown menu)
        text_special,
        debugging
    ]
        
    return all_page_content

def generate_dig_report_dataframe_combined(
        path_to_dig_results, 
    ):
    """
    Generates the dataframe that contains the data for the DIG report

    Parameters
    ==========
        path_to_dig_results: str,
            file path name to the combined dig report data
            
    Return
    ======
        pd.DataFrame
            dataframe for the combined dig report type
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

def gen_coding_region_app_component(
        dig_df,
        mutation_type,
        burden_type,
        p_val_type,
        display_bounds_value,
    ):  
    """
    
    """
    debugging = ""
    all_page_content = []
    dig_df = dig_df.copy()
    bad_mutation_burden_combination = [("indels_nonsynonymous_snvs", "sample_wise"), ("indels", "sample_wise")]
    
    if (mutation_type, burden_type) not in bad_mutation_burden_combination:
        # sort the table with respect to 'PVAL' column, smallest(1) -> largest(nth)
        dig_df = dig_df.sort_values(by='PVAL'+ "_"+ coding_region_mutation_type[mutation_type] + "_" + coding_region_burden_type[burden_type]) #+ "_" + scatterpoint_type[p_val_type])
        dig_df['RANK'] = np.array([i+1 for i in range(len(dig_df))])

    # checks if the user wants to display the bounds on the dig report plot
    if display_bounds_value:
        display_bounds_key = 'Yes'

    # defaults to not displaying lower/upper bounds on the dig report plot
    else:
        display_bounds_key = 'No'

    # MAKE SURE TO REDO THE CODING REGION VALUES, ADD IF STATEMENTS TO CHANGE WHAT GETS DISPLAYED BASED ON THE DROP DOWN MENU
    df_kept, volcano_fig, qq_fig, fig_mu, fig_sigma, dnds_fig, table_fig, text_special = generate_coding_region_report(dig_df, 
                                                                                                                       mutation_type, 
                                                                                                                       burden_type, 
                                                                                                                       display_bounds_key, 
                                                                                                                       p_val_type)
    dig_data_columns = []
        
    for clm_nm in DIG_CODING_REGION_REPORT_COLUMN_NAMES:
        dig_data_clm_dict = {}
        dig_data_clm_dict["name"] = clm_nm
        dig_data_clm_dict["id"] = clm_nm
        print("this is the column name: ", clm_nm)
        
        try: 
            if "PVAL" in clm_nm:
                dig_data_clm_dict["id"] = clm_nm + "_"+ coding_region_mutation_type[mutation_type] + "_" + coding_region_burden_type[burden_type] #+ "_" +scatterpoint_type[p_val_type]

            elif 'FDR' in clm_nm:
                dig_data_clm_dict["id"] = clm_nm + "_" + coding_region_mutation_type[mutation_type] + "_" + coding_region_burden_type[burden_type] + "_" + scatterpoint_type[p_val_type]
            
            elif 'LENGTH' in clm_nm:
                dig_data_clm_dict['id'] = 'GENE_' +  clm_nm

            # no modifications needed 
            elif 'dNdS_OBS' in clm_nm or 'dNdS_EXP' in clm_nm:
                dig_data_clm_dict["id"] = clm_nm 

            elif 'OBS' in clm_nm or 'EXP' in clm_nm: 
                dig_data_clm_dict["id"] = clm_nm + "_"+ coding_region_mutation_type[mutation_type]
            
            dig_data_columns.append(dig_data_clm_dict)

        # skip adding this column to the data dictionary 
        except KeyError as k:
            debugging = debugging + f"error: {k}/" + clm_nm + "_" + coding_region_mutation_type[mutation_type] + coding_region_burden_type[burden_type]
            print("ran into an issue with this column name: ", clm_nm + "_" + coding_region_mutation_type[mutation_type] + "_"+ coding_region_burden_type[burden_type])
            continue

    if len(df_kept):
        for column_dict in dig_data_columns:
            # gets the dig data column name
            column = column_dict["id"]

            # reformats the columns with float values in them
            if pd.api.types.is_float_dtype(df_kept[column]) and column != 'RANK':
                df_kept[column] = reformat_numbers(df_kept[column], format='{:.3E}')
            
        # ONLY GETTING THE FIRST 100 ROWS OF DATA TO DISPLAY IN THE TABLE
        # REMOVE THE DEBUGGING LATER!!!
        df_kept = df_kept[:100] 

    all_page_content = [
            df_kept.to_dict('records'),
            volcano_fig,
            qq_fig,
            fig_mu, 
            fig_sigma,
            dnds_fig,
            dig_data_columns,
            mutation_type,
            burden_type,
            p_val_type,

            # updates the dropdown options
            CODING_REGION_MUT_DROPDOWN,
            CODING_REGION_BUR_DROPDOWN,
            CODING_REGION_SCATTER_DROPDOWN, 

            # changing the plot size and the visibility of the plots
            SHOW_QQ_SMALL_PLOT_STYLE,   # displays the qq plot
            SHOW_VOLCANO_PLOT_STYLE,    # displays the volcano plot
            SHOW_THREE_PLOTS_STYLE,     # displays the fig mu plot
            SHOW_THREE_PLOTS_STYLE,     # displays the fig sigma plot
            SHOW_THREE_PLOTS_STYLE,     # displays the dnds plot

            text_special,
            debugging
        ]
    return all_page_content

def gen_non_coding_app_component(
    dig_df,
    mutation_type,
    burden_type,
    p_val_type,
    display_bounds_value,
):
    """
    
    """
    debugging=""
    all_page_content = []
    dig_df = dig_df.copy()
    bad_mutation_burden_combination = [("indels_snvs", "sample_wise"), ("indels", "sample_wise")]

    if (mutation_type, burden_type) not in bad_mutation_burden_combination:
        column_to_sort_by = ""
        
        if burden_type == "sample_wise":
            column_to_sort_by = 'PVAL'+ "_" + combined_mutation_type[mutation_type] + "_BURDEN"
        else:
            column_to_sort_by = 'PVAL'+ "_" + combined_mutation_type[mutation_type] + "_" + non_coding_region_burden_type[burden_type]
        dig_df = dig_df.sort_values(by=column_to_sort_by)
        dig_df['RANK'] = np.array([i+1 for i in range(len(dig_df))])

    # initialize the plots for the 3 prime str report
    qq_fig = go.Figure()
    fig_mu = go.Figure()
    fig_sigma = go.Figure()
    dnds_fig = go.Figure()
    volcano_fig = go.Figure()

    if display_bounds_value:
        display_bounds_key = 'Yes'
    
    else:
        display_bounds_key = 'No'

    df_kept, volcano_fig, qq_fig, fig_mu, fig_sigma, table_fig, text_special = generate_dig_non_coding_plots(dig_df, mutation_type, burden_type, display_bounds_key, p_val_type)
    dig_data_columns = []
    
    for clm_nm in DIG_NON_CODING_REPORT_COLUMN_NAMES:
        dig_data_clm_dict = {}
        dig_data_clm_dict["name"] = clm_nm
        dig_data_clm_dict["id"] = clm_nm

        try:
            if "PVAL" in clm_nm and burden_type != "sample_wise":
                # gets the column name corresponding to the mutation type and the scatterpoint type
                dig_data_clm_dict["id"] = clm_nm + "_" + combined_mutation_type[mutation_type] + "_" + non_coding_region_burden_type[burden_type]
            elif "PVAL" in clm_nm:
                dig_data_clm_dict["id"] = clm_nm + "_" + non_coding_region_burden_type[burden_type]

            elif 'EXP' in clm_nm:
                dig_data_clm_dict["id"] = clm_nm + "_"+ combined_mutation_type[mutation_type]

            elif 'SIZE' in clm_nm:
                dig_data_clm_dict['id'] = 'ELT_SIZE'

            elif 'FDR' in clm_nm:
                dig_data_clm_dict["id"] = clm_nm + "_" + combined_mutation_type[mutation_type] + "_" + non_coding_region_burden_type[burden_type] +  "_" + scatterpoint_type[p_val_type]
            
            elif 'LENGTH' in clm_nm:
                dig_data_clm_dict['id'] = 'GENE_' +  clm_nm

            elif 'OBS' in clm_nm and combined_mutation_type[mutation_type] != 'MUT' and combined_mutation_type[mutation_type] != 'dNdS':
                dig_data_clm_dict["id"] = clm_nm + "_"+ combined_mutation_type[mutation_type]

            elif 'OBS' in clm_nm and mutation_type == 'indels_snvs':
                dig_data_clm_dict["id"] = clm_nm + "_SAMPLES"    

            dig_data_columns.append(dig_data_clm_dict)
        except KeyError as k:
            print("experiencing key error with this column: ", clm_nm)
    
    if len(df_kept):
        for column_dict in dig_data_columns:
            # gets the dig data column name
            column = column_dict["id"]
            
            # reformats the columns with float values in them
            if pd.api.types.is_float_dtype(df_kept[column]) and column != 'RANK':
                df_kept[column] = reformat_numbers(df_kept[column], format='{:.3E}')

        # ONLY GETTING THE FIRST 100 ROWS OF DATA TO DISPLAY IN THE TABLE
        # REMOVE THE DEBUGGING LATER!!! 
        df_kept = df_kept[:100]

    all_page_content =  [
            df_kept.to_dict('records'),
            volcano_fig,
            qq_fig,
            fig_mu, 
            fig_sigma,
            dnds_fig,
            dig_data_columns,
            mutation_type,
            burden_type,
            p_val_type,

            # dropdown options
            COMBINED_MUT_DROPDOWN,
            CODING_REGION_BUR_DROPDOWN,
            CODING_REGION_SCATTER_DROPDOWN,

            # changing the plot size and the visibility of the plots
            SHOW_QQ_SMALL_PLOT_STYLE,   # displays the qq plot
            SHOW_VOLCANO_PLOT_STYLE,    # displays the volcano plot
            SHOW_TWO_PLOTS_STYLE,       # displays the fig mu plot
            SHOW_TWO_PLOTS_STYLE,       # displays the fig sigma plot
            HIDE_PLOT_STYLE,            # hides the dnds plot

            text_special,
            debugging
        ]
    return all_page_content