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

from SNVReviewers.AppComponents.utils import reformat_numbers

RESULTS_DROPDOWN_OPTIONS = [
    {'label':'All', 'value': 'all'},
    {'label': 'First 30', 'value': 'First 30'},
    {'label': 'First 20', 'value': 'First 20'},
    {'label': 'First 10', 'value': 'First 10'},
    {'label': 'None', 'value': 'None'}
]

COMPARISON_DROPDOWN_OPTIONS = [
    {'label': "MutSig2 vs dNdScv", 'value':'MutSig2CV vs dNdScv'}, 
    {'label': "MutSig2 vs DIG", 'value':"MutSig2CV vs DIG"}, 
    {'label': "dNdScv vs DIG", 'value':"dNdScv vs DIG"}
]

SUMMARY_DROPDOWN_OPTIONS = [
    {'label': "MutSig2, dNdScv, and DIG", 'value':"MutSig2CV, dNdScv, and DIG"}, 
    {'label': "MutSig2 and dNdScv only", 'value':"MutSig2CV and dNdScv only"}, 
    {'label': "MutSig2 and DIG only", 'value':"MutSig2CV and DIG only"},
    {'label': "dNdScv and DIG only", 'value':"dNdScv and DIG only"}, 
    {'label': "MutSig2 only", 'value':"MutSig2CV only"}, 
    {'label': "dNdScv only", 'value':"dNdScv only"}, 
    {'label': "DIG only", 'value':"DIG only"}
]

QQ_PLOT_RESULT_STYLE = {"width":"1200px"}
DNDS_GLOBAL_RESULT_STYLE = {'display':'block', 'width':'600px'}
DNDS_MIS_RESULT_STYLE = {'display':'block', 'width':'600px'}
DNDS_TRUNC_RESULT_STYLE = {'display':'block', 'width':'600px'}
DNDS_COMPARISON_STYLE = {"width":"1200px"}
HIDE_PLOTS_STYLE = {"display":"none"}

HIDE_TABLE_STYLE = {"display":"none"}
SHOW_TABLE_STYLE = {
                    'width': '100%',  # Make the table width responsive
                    'maxWidth': '100%',  # Ensure it doesn’t go beyond the screen width
                    'overflowX': 'auto',  # Allow horizontal scroll if necessary
                    }

DNDSCV_REPORT_COLUMN_NAMES = ["RANK", "GENE", "N_SYN", "N_MIS", "N_NON", "N_SPL", "N_IND",
                              "dNdS_MIS", "dNdS_NON", "dNdS_SPL", "dNdS_IND", "PVAL_MIS",
                              "PVAL_TRUNC", "PVAL_IND", "PVAL", "FDR", "CGC", "PANCAN"]

DNDSCV_SUMMARY_COLUMN_NAMES = ['RANK', 'GENE', 'SIZE_coding', 'PVAL_MutSig2', 'PVAL_dNdScv', 'PVAL_DIG', 
                               'PVAL_comb', 'FDR_MutSig2', 'FDR_dNdScv', 'FDR_DIG', 'FDR_comb', 'SIG_MutSig2',
                               'SIG_dNdScv', 'SIG_DIG', 'CGC', 'PANCAN']

DNDS_RESULTS_DROPDOWM_LABEL = 'Select Number of Significant Gene Labels to Display:'
DNDS_COMPARISON_DROPDOWM_LABEL = 'Tools to compare:'
DNDS_SUMMARY_DROPDOWM_LABEL = 'List genes significant with:'

# dNdScv dataframe and plot generation
import argparse
import pandas as pd
import plotly.graph_objects as go
import numpy as np
import json
from SNVReviewers.AppComponents.utils import generate_dnds_report, gen_dnds_summary_table, gen_dnds_comparison_plot

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
    summary_table = go.Figure()
    comparison_table = go.Figure()

    # for clm in dnd_df_plot:
    #     debugging = debugging + " " + clm

    qq_fig, fig_dnds_global, fig_dnds_mis, fig_dnds_tru, df_plot = generate_dnds_report(dnd_df_plot, dnd_df_merged,dnd_df_global, num_gene_values)

    for clm_nm in DNDSCV_REPORT_COLUMN_NAMES:
        
        # reformats the columns with float values in them
        if pd.api.types.is_float_dtype(dnd_df_plot[clm_nm]):
            dnd_df_plot[clm_nm] = reformat_numbers(dnd_df_plot[clm_nm], format='{:.3E}')
    
    # NEED TO REFORMAT THE dnd_df_plot dataframe for 2 point decimal precision
    all_page_content = [
            # data to be displayed in the tables for the results dndscv report
            dnd_df_plot.to_dict('records'),
            comparison_table,
            comparison_table,
            comparison_table,
            summary_table, 

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
    dnd_df_plot,
    dnd_df_merged,
    dnd_df_global,
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
    empty_summary_table = go.Figure()

    comparison_fig, [comaprison_table1, comaprison_table2, comaprison_table3] = gen_dnds_comparison_plot(dnd_df_comparison, dropdown_menu_value)

    all_page_content = [
        # data to be displayed in the tables for the comparison dndscv report
        empty_result_table,
        comaprison_table1,
        comaprison_table2,
        comaprison_table3,
        empty_summary_table,

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

        # comparison dndscv dropdown menu values and labels
        dnd_radio_item_selection,
        COMPARISON_DROPDOWN_OPTIONS,
        DNDS_COMPARISON_DROPDOWM_LABEL,
        dropdown_menu_value,

        debugging
    ]

    return all_page_content

def gen_dndscv_summary_app_component(
    dnd_df_plot,
    dnd_df_merged,
    dnd_df_global,
    dnd_df_comparison,
    dropdown_menu_value,
    dnd_radio_item_selection       
):
    """
    """
    debugging= ""
    all_page_content = []
    empty_figure = go.Figure()
    
    summary_table1 = pd.DataFrame().to_dict('records')

    # FINISH THIS LATER!!!
    print("before I try to generate dnds summary table!!")
    # summary_table_dict, all_titles = gen_dnds_summary_table(dnd_df_comparison, dropdown_menu_value)
    # summary_table2 = gen_dnds_summary_table(dnd_df_comparison, dropdown_menu_value)

    summary_tables, all_titles = gen_dnds_summary_table(dnd_df_comparison, dropdown_menu_value)

    for title in all_titles:
        debugging = debugging +"/" + title
    summary_table_combined = go.Figure()
    # summary_table_combined = summary_tables["MutSig2CV, dNdScv, and DIG"]
    print("this is dropdown_menu_value: ", dropdown_menu_value)
    # print("these are the keys in the summary table dict: ", list(summary_table_dict.keys()))
    if dropdown_menu_value != "MutSig2CV, dNdScv, and DIG":

        summary_table2 = summary_tables[dropdown_menu_value]
    else:
        summary_table2 = summary_tables["dNdScv and DIG only"]

    comparison_table = go.Figure()

    # table5 = dnd_df_merged.to_dict('records'),

    all_page_content = [
        # tables to display for the summary dNdScv report
        summary_table1,
        summary_table_combined,
        comparison_table,
        comparison_table,
        summary_table2,

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

        # summary dndscv dropdown menu values and labels
        dnd_radio_item_selection,
        SUMMARY_DROPDOWN_OPTIONS,
        DNDS_SUMMARY_DROPDOWM_LABEL,
        dropdown_menu_value,

        debugging
    ]

    return all_page_content