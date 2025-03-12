import plotly.graph_objects as go
import numpy as np
import pandas as pd

from SNVReviewers.AppComponents.utils import generate_combined_dig_report_plots, combined_mutation_type, combined_burden_plot_type, scatterpoint_type
from SNVReviewers.AppComponents.utils import generate_coding_region_report, coding_region_mutation_type, coding_region_burden_type
from SNVReviewers.AppComponents.utils import generate_dig_non_coding_plots, non_coding_region_burden_type, reformat_numbers

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

def gen_combined_app_component(
    dig_df,
    mutation_type,
    burden_type,
    p_val_type,
    display_bounds_value,
):
    """
    
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