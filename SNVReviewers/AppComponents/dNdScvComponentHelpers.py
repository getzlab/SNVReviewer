
import plotly.graph_objects as go
import pandas as pd
import numpy as np
from SNVReviewers.AppComponents.utils import reformat_numbers
import pandas as pd
import numpy as np
from SNVReviewers.AppComponents.utils import generate_dnds_report, gen_dnds_summary_table, gen_dnds_comparison_plot, gen_table_conditional_styling

RESULTS_DROPDOWN_OPTIONS = [
    {'label':'All', 'value': 'all'},
    {'label': 'First 30', 'value': 'First 30'},
    {'label': 'First 20', 'value': 'First 20'},
    {'label': 'First 10', 'value': 'First 10'},
    {'label': 'None', 'value': 'None'}
]

COMPARISON_DROPDOWN_OPTIONS = [
    {'label': "MutSig2CV vs dNdScv", 'value':'MutSig2CV vs dNdScv'}, 
    {'label': "MutSig2CV vs DIG", 'value':"MutSig2CV vs DIG"}, 
    {'label': "dNdScv vs DIG", 'value':"dNdScv vs DIG"}
]

SUMMARY_DROPDOWN_OPTIONS = [
    {'label': "MutSig2CV, dNdScv and DIG", 'value':"MutSig2CV, dNdScv and DIG"}, 
    {'label': "MutSig2CV and dNdScv only", 'value':"MutSig2CV and dNdScv only"}, 
    {'label': "MutSig2CV and DIG only", 'value':"MutSig2CV and DIG only"},
    {'label': "dNdScv and DIG only", 'value':"dNdScv and DIG only"}, 
    {'label': "MutSig2CV only", 'value':"MutSig2CV only"}, 
    {'label': "dNdScv only", 'value':"dNdScv only"}, 
    {'label': "DIG only", 'value':"DIG only"}
]

# different plot styles
QQ_PLOT_RESULT_STYLE = {"width":"1200px"}
DNDS_GLOBAL_RESULT_STYLE = {'display':'block', 'width':'600px'}
DNDS_MIS_RESULT_STYLE = {'display':'block', 'width':'600px'}
DNDS_TRUNC_RESULT_STYLE = {'display':'block', 'width':'600px'}
DNDS_COMPARISON_STYLE = {"width":"1200px"}
HIDE_PLOTS_STYLE = {"display":"none"}

# different table styles
HIDE_TABLE_STYLE = {"display":"none"}
SHOW_TABLE_STYLE = {
                    'width': '100%',  # Make the table width responsive
                    'maxWidth': '100%',  # Ensure it doesn’t go beyond the screen width
                    'overflowX': 'auto',  # Allow horizontal scroll if necessary
                    }

# column names to be displayed on the dash table depending on report type
DNDSCV_REPORT_COLUMN_NAMES = ["RANK", "GENE", "N_SYN", "N_MIS", "N_NON", "N_SPL", "N_IND",
                              "dNdS_MIS", "dNdS_NON", "dNdS_SPL", "dNdS_IND", "PVAL_MIS",
                              "PVAL_TRUNC", "PVAL_IND", "PVAL", "FDR", "CGC", "PANCAN"]

DNDSCV_COMPARISON_COLUMN_NAMES = ['RANK', 'GENE', 'SIZE_coding', 'PVAL_MutSig2', 'PVAL_dNdScv', 'PVAL_DIG', 
                               'PVAL_comb', 'FDR_MutSig2', 'FDR_dNdScv', 'FDR_DIG', 'FDR_comb', 'SIG_MutSig2',
                               'SIG_dNdScv', 'SIG_DIG', 'CGC', 'PANCAN']

DNDSCV_SUMMARY_COLUMN_NAMES = ['RANK', 'GENE', 'SIZE_coding', 'PVAL_MutSig2', 'PVAL_dNdScv', 'PVAL_DIG', 
                               'PVAL_comb', 'FDR_MutSig2', 'FDR_dNdScv', 'FDR_DIG', 'FDR_comb', 'SIG_MutSig2',
                               'SIG_dNdScv', 'SIG_DIG', 'CGC', 'PANCAN']

# dropdown label for the different report types
DNDS_RESULTS_DROPDOWM_LABEL = 'Select Number of Significant Gene Labels to Display:'
DNDS_COMPARISON_DROPDOWM_LABEL = 'Tools to compare:'
DNDS_SUMMARY_DROPDOWM_LABEL = 'List genes significant with:'
NO_CONDITIONAL_STYLING = []

# dNdScv dataframe and plot generation
def gen_dndscv_results_app_component(
    dnd_df_plot,
    dnd_df_merged,
    dnd_df_global,
    num_gene_values,
    dnd_radio_item_selection
):
    """ 
    """
    debugging= ""
    all_page_content = []
    dnd_df_plot = dnd_df_plot.copy()
    dnd_df_merged = dnd_df_merged.copy()
    dnd_df_global = dnd_df_global.copy()
    summary_table = pd.DataFrame().to_dict("records")
    comparison_table = pd.DataFrame().to_dict("records")

    qq_fig, fig_dnds_global, fig_dnds_mis, fig_dnds_tru, df_plot, signficant_boolean_values = generate_dnds_report(dnd_df_plot, dnd_df_merged,dnd_df_global, num_gene_values)
    signficant_idxs = np.where(signficant_boolean_values)[0]

    for clm_nm in DNDSCV_REPORT_COLUMN_NAMES:
        
        # reformats the columns with float values in them
        if pd.api.types.is_float_dtype(dnd_df_plot[clm_nm]):
            dnd_df_plot[clm_nm] = reformat_numbers(dnd_df_plot[clm_nm], format='{:.3E}')

    dnds_comparison_column_name = [
        {"name":i, 
         "id":i} for i in DNDSCV_COMPARISON_COLUMN_NAMES
        ]
    dnds_summary_column_names = [
        {"name":i, 
         "id":i} for i in DNDSCV_SUMMARY_COLUMN_NAMES
    ]

    conditional_styling = gen_table_conditional_styling(signficant_idxs)
    
    all_page_content = [
            # data to be displayed in the tables for the results dndscv report
            dnd_df_plot.to_dict('records'),
            comparison_table,
            comparison_table,
            comparison_table,
            summary_table, 

            # dnds comparison and summary report type column names
            dnds_comparison_column_name,
            dnds_comparison_column_name,
            dnds_comparison_column_name,
            dnds_summary_column_names,

            # result dndscv report figures
            qq_fig,
            fig_dnds_global,
            fig_dnds_mis,
            fig_dnds_tru,
            go.Figure(), # do not display comparison plot

            # result dndscv report figure style
            QQ_PLOT_RESULT_STYLE,
            DNDS_GLOBAL_RESULT_STYLE,
            DNDS_MIS_RESULT_STYLE,
            DNDS_TRUNC_RESULT_STYLE, 
            HIDE_PLOTS_STYLE,

            # result dndscv report table style
            SHOW_TABLE_STYLE,
            HIDE_TABLE_STYLE,
            HIDE_TABLE_STYLE,
            HIDE_TABLE_STYLE,
            HIDE_TABLE_STYLE,

            # result dndscv report table conditional styling
            conditional_styling,
            NO_CONDITIONAL_STYLING,
            NO_CONDITIONAL_STYLING,
            NO_CONDITIONAL_STYLING,
            NO_CONDITIONAL_STYLING,
            
            # result dndscv dropdown menu values and labels
            dnd_radio_item_selection,
            RESULTS_DROPDOWN_OPTIONS,
            DNDS_RESULTS_DROPDOWM_LABEL,
            num_gene_values,

            debugging
        ]

    # NEED TO ASK DAVID WHICH COLUMN TO SORT THE dnd_df by
    # dnd_df = dnd_df.sort_values(by='PVAL'+ "_"+ combined_mutation_type[mutation_type] + "_" + scatterpoint_type[p_val_type])
    # dnd_df['RANK'] = np.array([i+1 for i in range(len(dnd_df))])

    return all_page_content

def gen_dndscv_comparison_app_component(
    dnd_df_comparison,
    dropdown_menu_value,
    dnd_radio_item_selection           
):
    """ 
    """
    debugging= ""
    all_page_content = []
    empty_figure = go.Figure()
    comparison_fig = go.Figure()
    empty_result_table = pd.DataFrame().to_dict('records')
    empty_summary_table = pd.DataFrame().to_dict("records")

    comparison_fig, [df_plot1, df_plot2, df_plot3], df_significance = gen_dnds_comparison_plot(dnd_df_comparison, dropdown_menu_value)

    dnds_comparison_column_names1 =[
                            {"name": i,
                            "id": i} for i in list(df_plot1.columns)
                        ]
    dnds_comparison_column_names2 =[
                            {"name": i,
                            "id": i} for i in list(df_plot2.columns)
                        ]
    dnds_comparison_column_names3 =[
                            {"name": i,
                            "id": i} for i in list(df_plot3.columns)
                        ]
    
    dnds_summary_column_names = [
        {"name":i, 
         "id":i} for i in DNDSCV_SUMMARY_COLUMN_NAMES
    ]

    for clm_nm in df_plot1.columns:
        debugging = debugging + "/" + clm_nm

    significant_idxs = np.where(((df_significance['Significant'].str.contains('both')) | (df_significance['Significant'].str.contains('only'))) )[0]
    # significant_idxs2 = np.where(~df_plot2['Significant'].str.contains('neither'))[0]
    # significant_idxs3 = np.where(~df_plot3['Significant'].str.contains('neither'))[0]
    debugging = f"{type(significant_idxs)}"
    print("")
    
    comparison_table_conditional_styling1 = gen_table_conditional_styling(significant_idxs)
    comparison_table_conditional_styling2 = gen_table_conditional_styling(significant_idxs)
    comparison_table_conditional_styling3 = gen_table_conditional_styling(significant_idxs)

    # DEBUGGING MAKING SURE THAT THE HIGHLIGHTING/BOLDING WORKS!!!
    comparison_table_conditional_styling1 = [
            # bolds rows corresponding to genes that are significant
            {
                'if': {
                    'row_index': [0, 1, 2, 3]
                },
                'fontWeight': 'bold',  # Bold the font
            },

            # highlights rows corresponding to genes that are significant
            {
                'if': {
                    'row_index': [0, 1, 2, 3]
                },
                'backgroundColor': 'yellow',  
            },
        ]

    all_page_content = [
        # data to be displayed in the tables for the comparison dndscv report
        empty_result_table,
        df_plot1.to_dict("records"),
        df_plot2.to_dict("records"),
        df_plot3.to_dict("records"),
        empty_summary_table,

        # dnds comparison and summary report type column names
        dnds_comparison_column_names1,
        dnds_comparison_column_names2,
        dnds_comparison_column_names3,
        dnds_summary_column_names,

        # comparison dndscv report figures
        empty_figure,
        empty_figure, 
        empty_figure,
        empty_figure,
        comparison_fig,
        
        # comparison dndscv report figure style
        HIDE_PLOTS_STYLE, 
        HIDE_PLOTS_STYLE,
        HIDE_PLOTS_STYLE,
        HIDE_PLOTS_STYLE,
        DNDS_COMPARISON_STYLE,

        # comparison dndscv report table style
        HIDE_TABLE_STYLE,
        SHOW_TABLE_STYLE,
        SHOW_TABLE_STYLE,
        SHOW_TABLE_STYLE,
        HIDE_TABLE_STYLE,

        # comparison dndscv report table conditional styling
        NO_CONDITIONAL_STYLING,
        comparison_table_conditional_styling1,
        comparison_table_conditional_styling2,
        comparison_table_conditional_styling3,
        NO_CONDITIONAL_STYLING,

        # comparison dndscv dropdown menu values and labels
        dnd_radio_item_selection,
        COMPARISON_DROPDOWN_OPTIONS,
        DNDS_COMPARISON_DROPDOWM_LABEL,
        dropdown_menu_value,

        debugging
    ]

    return all_page_content

def gen_dndscv_summary_app_component(
    dnd_df_comparison,
    dropdown_menu_value,
    dnd_radio_item_selection       
):
    """
    """
    debugging= ""
    all_page_content = []
    empty_figure = go.Figure()
    empty_comparison_table = pd.DataFrame().to_dict("records")
    summary_table1 = pd.DataFrame().to_dict('records')

    df_summary_tables = gen_dnds_summary_table(dnd_df_comparison, dropdown_menu_value)
    summary_table_combined = df_summary_tables["MutSig2CV, dNdScv and DIG"]

    print("this is dropdown_menu_value: ", dropdown_menu_value)

    summary_table2 = df_summary_tables[dropdown_menu_value]

    dnds_comparison_column_names = [
        {"name":i, 
         "id":i} for i in DNDSCV_COMPARISON_COLUMN_NAMES
    ]

    dnds_summary_column_names = [
        {"name":i, 
         "id":i} for i in list(summary_table_combined.columns)
    ]
    
    for clm_nm in summary_table_combined.columns:
        if 'FDR' in clm_nm:
            summary_table_conditional_styling = gen_table_conditional_styling(clm_nm, 0.1)
     

    all_page_content = [
        # tables to display for the summary dNdScv report
        summary_table1,
        summary_table_combined.to_dict("records"),
        empty_comparison_table,
        empty_comparison_table,
        summary_table2.to_dict("records"),

        # dnds comparison and summary report type column names
        dnds_summary_column_names,
        dnds_comparison_column_names,
        dnds_comparison_column_names,
        dnds_summary_column_names,

        # no figures to display for the summary dNdScv report
        empty_figure,
        empty_figure, 
        empty_figure,
        empty_figure,
        empty_figure,

        # summary dndscv report figure style
        HIDE_PLOTS_STYLE,
        HIDE_PLOTS_STYLE,
        HIDE_PLOTS_STYLE,
        HIDE_PLOTS_STYLE,
        HIDE_PLOTS_STYLE,

        # summary dndscv report table style
        HIDE_TABLE_STYLE,
        SHOW_TABLE_STYLE,
        HIDE_TABLE_STYLE,
        HIDE_TABLE_STYLE,
        SHOW_TABLE_STYLE,

        # summary dndscv report table conditional styling
        NO_CONDITIONAL_STYLING,
        summary_table_conditional_styling,
        NO_CONDITIONAL_STYLING,
        NO_CONDITIONAL_STYLING,
        summary_table_conditional_styling,

        # summary dndscv dropdown menu values and labels
        dnd_radio_item_selection,
        SUMMARY_DROPDOWN_OPTIONS,
        DNDS_SUMMARY_DROPDOWM_LABEL,
        dropdown_menu_value,

        debugging
    ]

    return all_page_content