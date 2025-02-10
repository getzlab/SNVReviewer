from dash import dcc, html
from dash.dependencies import Input, Output
import dash_bootstrap_components as dbc
import plotly.graph_objects as go
import dash_daq as daq

from AnnoMate.ReviewDataApp import AppComponent
from AnnoMate.DataTypes.GenericData import GenericData

from dash import dcc, html, dash_table
from dash.dependencies import Input, Output
import dash_bootstrap_components as dbc
import plotly.graph_objects as go

from AnnoMate.Data import Data, DataAnnotation
from AnnoMate.ReviewDataApp import ReviewDataApp, AppComponent
from AnnoMate.DataTypes.GenericData import GenericData
from cnv_suite.visualize import plot_acr_interactive

from rpy2.robjects import r, pandas2ri
import rpy2.robjects as robjects
import os
import pickle
from typing import Union, List, Dict
import sys
from cnv_suite import calc_cn_levels
import pandas as pd
import numpy as np

# from SNVReviewers.AppComponents.DIGAppComponent_CodingRegion import gen_dig_app_coding_region_component_layout, DIG_CODING_REGION_REPORT_COLUMN_NAMES
from SNVReviewers.AppComponents.utils import generate_combined_dig_report_plots, combined_mutation_type, scatterpoint_type
from SNVReviewers.AppComponents.utils import generate_coding_region_report, coding_region_mutation_type, coding_region_burden_type
from SNVReviewers.AppComponents.utils import generate_dig_non_coding_plots

DIG_CODING_REGION_REPORT_COLUMN_NAMES = ["RANK", "GENE", 'CHROM', 'LENGTH', "FDR", "PVAL", "OBS", 
                                         'EXP', 'MU', 'SIGMA', 'dNdS_OBS', 'dNdS_EXP', 'FLAG', 'CGC', 'PANCAN']

DIG_REPORT_COMBINED_COLUMN_NAMES = ["RANK", "GENE", "FDR", "PVAL", "PVAL_coding", "PVAL_promoter", "PVAL_5utr", 
                           "SIZE_coding", "SIZE_promoter", "SIZE_5utr", "SIZE_3utr", "CGC", "PANCAN"]

DIG_NON_CODING_REPORT_COLUMN_NAMES = ["RANK", "GENE", 'SIZE', "PVAL", 
                                    #   "FDR",  # Figure out why this column is not being displayed!!
                                      "OBS", 
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

DIG_REPORT_VALUES = ["Combined", "Coding region", "Promoter region", "5-prime UTRs", "3-prime UTRs"]
DIG_DATAFRAME_IDX = 0
DND_DATAFRAME_IDX = 1
MUTSIG_DATAFRAME_IDX = 2
SNV_DATA_COLUMN_NAME = "snv_data"
DIG_LABEL_IDX = 1

def gen_combined_app_component(
    dig_df,
    dig_type_selection, # radio item selection
    mutation_type,
    burden_type,
    p_val_type,
    display_toggle_value,
    display_label_value
):
    """
    
    """
    all_page_content = []
    # dig_df = data.df[SNV_DATA_COLUMN_NAME][0][DIG_DATAFRAME_IDX].copy() # gets the dig report data for a specific cohort
    dig_df = dig_df.copy()
    dig_df = dig_df.sort_values(by='PVAL'+ "_"+ combined_mutation_type[mutation_type] + "_" + scatterpoint_type[p_val_type])
    dig_df['RANK'] = np.array([i+1 for i in range(len(dig_df))])

    display_bounds = display_toggle_value # whether to display the bounds
    
    qq_fig = go.Figure()
    fig_mu = go.Figure()
    fig_sigma = go.Figure()
    dnds_fig = go.Figure()
    volcano_fig = go.Figure()

    # checks if the user wants to display the bounds on the dig report plot
    if display_toggle_value:
        display_bounds = 'Yes'

    # defaults to not displaying lower/upper bounds on the dig report plot
    else:
        display_bounds = 'No'
    
    if display_label_value:
        display_labels_key = 'Yes'
    
    else:
        display_labels_key = 'No'

    qq_fig, table_fig, text_special = generate_combined_dig_report_plots(dig_df, mutation_type, 
                                                                            burden_type, display_bounds, 
                                                                            display_labels_key, p_val_type)
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

        # skip the rank column
        if column == 'RANK':
            continue

        format='{:.3E}'

        # rounds all the values in the FDR and PVAL columns to 4 significant digits
        if 'FDR' in column or 'PVAL' in column:
            dig_df[column] = [format.format(value) for value in dig_df[column]]
        
    # ONLY GETTING THE FIRST 100 ROWS OF DATA TO DISPLAY IN THE TABLE
    # REMOVE THE DEBUGGING LATER!!!   
    dig_df = dig_df[:100]

    all_page_content = [
                dig_df.to_dict('records'),
                dig_type_selection,
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
                {'width':'1200px'},
                {'display':'none'}, # hides the volcano plot
                {'display':'none'}, # hides the volcano plot
                {'display':'none'}, # hides the volcano plot
                {'display':'none'}, # hides the volcano plot
                # Warning text message (only generated if you get a specific mutation and burden type combination on dropdown menu)
                text_special,
            ]
        
    return all_page_content

def gen_coding_region_app_component(
        dig_df,
        dig_type_selection, # radio item selection
        mutation_type,
        burden_type,
        p_val_type,
        display_toggle_value,
        display_label_value
    ):  
    """
    
    """
    all_page_content = []
    # dig_df = data.df[SNV_DATA_COLUMN_NAME][0][DIG_DATAFRAME_IDX].copy() # gets the dig report data for a specific cohort
    dig_df = dig_df.copy()
    # sort the table with respect to 'PVAL' column, smallest(1) -> largest(nth)
    dig_df = dig_df.sort_values(by='PVAL'+ "_"+ coding_region_mutation_type[mutation_type] + "_" + coding_region_burden_type[burden_type]) #+ "_" + scatterpoint_type[p_val_type])
    dig_df['RANK'] = np.array([i+1 for i in range(len(dig_df))])

    # debugging_component = ""    
    qq_fig = go.Figure()
    
    if display_label_value:
        display_labels_key = 'Yes'
    
    else:
        display_labels_key = 'No'

    # MAKE SURE TO REDO THE CODING REGION VALUES, ADD IF STATEMENTS TO CHANGE WHAT GETS DISPLAYED BASED ON THE DROP DOWN MENU
    df_kept, volcano_fig, qq_fig, fig_mu, fig_sigma, dnds_fig, table_fig, text_special = generate_coding_region_report(dig_df, mutation_type, burden_type, display_labels_key, p_val_type)
    df_kept_new = df_kept[['dNdS_OBS', 'dNdS_EXP', 'GENE']]
    dig_data_columns = []
    
    for clm_nm in DIG_CODING_REGION_REPORT_COLUMN_NAMES:
        dig_data_clm_dict = {}
        dig_data_clm_dict["name"] = clm_nm
        dig_data_clm_dict["id"] = clm_nm

        if "PVAL" in clm_nm:
            # gets the column name corresponding to the mutation type and the scatterpoint type
            dig_data_clm_dict["id"] = clm_nm + "_"+ coding_region_mutation_type[mutation_type] + "_" + coding_region_burden_type[burden_type] #+ "_" +scatterpoint_type[p_val_type]

        elif 'FDR' in clm_nm:
            dig_data_clm_dict["id"] = clm_nm + "_"+ coding_region_mutation_type[mutation_type] + "_" + scatterpoint_type[p_val_type]
        
        elif 'LENGTH' in clm_nm:
            dig_data_clm_dict['id'] = 'GENE_' +  clm_nm

        elif 'OBS' in clm_nm and coding_region_mutation_type[mutation_type] != 'MUT' and coding_region_mutation_type[mutation_type] != 'dNdS':
            dig_data_clm_dict["id"] = clm_nm + "_"+ coding_region_mutation_type[mutation_type]

        dig_data_columns.append(dig_data_clm_dict)
    
    for column_dict in dig_data_columns:
        # gets the dig data column name
        column = column_dict["id"]

        # skip the rank column
        if column == 'RANK':
            continue

        format='{:.3E}'

        # rounds all the values in the FDR and PVAL columns to 4 significant digits
        if 'FDR' in column or 'PVAL' in column or 'MU' in column or 'SIGMA' in column: #or 'dNdS_OBS' in column or 'dNdS_EXP' in column:
            dig_df[column] = [format.format(value) for value in dig_df[column]]
        
    # get the coding region working plots working!!
    # get the display bounds selection tool working 

    # ONLY GETTING THE FIRST 100 ROWS OF DATA TO DISPLAY IN THE TABLE
    # REMOVE THE DEBUGGING LATER!!!   
    dig_df = dig_df[:100]

    all_page_content = [
            dig_df.to_dict('records'),
            dig_type_selection,
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
            {'width':'600px'},
            {'display':'inline-block', 'width':'600px'}, # displays the volcano plot
            {'display':'block', 'width':'350px'}, # displays the fig mu plot
            {'display':'block', 'width':'350px'},  # displays the fig sigma plot
            {'display':'block', 'width':'350px'},  # displays the dnds plot

            text_special,
        ]
    return all_page_content

def gen_non_coding_app_component(
    dig_df,
    dig_type_selection,
    mutation_type,
    burden_type,
    p_val_type,
    display_toggle_value,
    display_label_value
):
    """
    
    """
    all_page_content = []
    # CHANGE THE INPUT, SO INSTEAD OF data, MAKE INPUT dig_df
    dig_df = dig_df.copy()
    # dig_df = data.df[SNV_DATA_COLUMN_NAME][0][DIG_DATAFRAME_IDX].copy() # gets the dig report data for a specific cohort
    dig_df = dig_df.sort_values(by='PVAL'+ "_"+ combined_mutation_type[mutation_type] + "_" + coding_region_burden_type[burden_type])
    dig_df['RANK'] = np.array([i+1 for i in range(len(dig_df))])

    # initialize the plots for the 3 prime str report
    qq_fig = go.Figure()
    fig_mu = go.Figure()
    fig_sigma = go.Figure()
    dnds_fig = go.Figure()
    volcano_fig = go.Figure()

    if display_label_value:
        display_labels_key = 'Yes'
    
    else:
        display_labels_key = 'No'

    df_kept, volcano_fig, qq_fig, fig_mu, fig_sigma, table_fig, text_special = generate_dig_non_coding_plots(dig_df, mutation_type, burden_type, display_labels_key, p_val_type)
    # df_kept = df_kept.copy()
    dig_data_columns = []
    
    for clm_nm in DIG_NON_CODING_REPORT_COLUMN_NAMES:
        dig_data_clm_dict = {}
        dig_data_clm_dict["name"] = clm_nm
        dig_data_clm_dict["id"] = clm_nm

        if "PVAL" in clm_nm:
            # gets the column name corresponding to the mutation type and the scatterpoint type
            dig_data_clm_dict["id"] = clm_nm + "_"+ combined_mutation_type[mutation_type] + "_" + coding_region_burden_type[burden_type]

        elif 'EXP' in clm_nm:
            dig_data_clm_dict["id"] = clm_nm + "_"+ combined_mutation_type[mutation_type]

        elif 'SIZE' in clm_nm:
            dig_data_clm_dict['id'] = 'ELT_SIZE'

        # elif 'FDR' in clm_nm:
        #     dig_data_clm_dict["id"] = clm_nm + "_"+ combined_mutation_type[mutation_type] + "_" + scatterpoint_type[p_val_type]
        
        elif 'LENGTH' in clm_nm:
            dig_data_clm_dict['id'] = 'GENE_' +  clm_nm

        elif 'OBS' in clm_nm and combined_mutation_type[mutation_type] != 'MUT' and combined_mutation_type[mutation_type] != 'dNdS':
            dig_data_clm_dict["id"] = clm_nm + "_"+ combined_mutation_type[mutation_type]

        dig_data_columns.append(dig_data_clm_dict)
    
    for column_dict in dig_data_columns:
        # gets the dig data column name
        column = column_dict["id"]

        # skip the rank column
        if column == 'RANK':
            continue

        format='{:.3E}'

        # rounds all the values in the FDR and PVAL columns to 4 significant digits
        if 'FDR' in column or 'PVAL' in column or 'MU' in column or 'SIGMA' in column:
            dig_df[column] = [format.format(value) for value in dig_df[column]]

    # get the coding region working plots working!!
    # get the display bounds selection tool working 

    # ONLY GETTING THE FIRST 100 ROWS OF DATA TO DISPLAY IN THE TABLE
    # REMOVE THE DEBUGGING LATER!!! 
    dig_df = dig_df[:100]
    # dig_kept = dig_kept[:100]

    all_page_content =  [
            dig_df.to_dict('records'),
            # df_kept.to_dict('records'),
            dig_type_selection,
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
            {'width':'600px'},
            {'display':'inline-block', 'width':'600px'}, # displays the volcano plot
            {'display':'block', 'width':'600px'}, # displays the fig mu plot
            {'display':'block', 'width':'600px'},  # displays the fig sigma plot
            {'display':'none'},  # displays the dnds plot

            text_special,
        ]
    return all_page_content