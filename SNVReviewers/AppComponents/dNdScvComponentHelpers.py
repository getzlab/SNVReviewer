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
    {'label': "MutSig2, dNdScv, and DIG", 'value':"MutSig2, dNdScv, and DIG"}, 
    {'label': "MutSig2 and dNdScv only", 'value':"MutSig2 and dNdScv only"}, 
    {'label': "MutSig2 and DIG only", 'value':"MutSig2 and DIG only"},
    {'label': "dNdScv and DIG only", 'value':"dNdScv and DIG only"}, 
    {'label': "MutSig2 only", 'value':"MutSig2 only"}, 
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

DNDS_RESULTS_DROPDOWM_LABEL = 'Select Number of Significant Gene Labels to Display:'
DNDS_COMPARISON_DROPDOWM_LABEL = 'Tools to compare:'
DNDS_SUMMARY_DROPDOWM_LABEL = 'List genes significant with:'

# dNdScv dataframe and plot generation
import argparse
import pandas as pd
import plotly.graph_objects as go
import numpy as np
import json
from SNVReviewers.AppComponents.utils import generate_dnds_report

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
    table = pd.DataFrame().to_dict('records')
    # table3 = pd.DataFrame().to_dict('records')


    for clm in dnd_df_plot:
        debugging = debugging + " " + clm
    # figure1 = go.Figure()

    qq_fig, fig_dnds_global, fig_dnds_mis, fig_dnds_tru, df_plot = generate_dnds_report(dnd_df_plot, dnd_df_merged,dnd_df_global, num_gene_values)

    # NEED TO REFORMAT THE dnd_df_plot dataframe for 2 point decimal precision
    all_page_content = [
            dnd_df_plot.to_dict('records'),
            table,
            table,
            table,
            table, 

            # dndscv report figures
            qq_fig,
            fig_dnds_global,
            fig_dnds_mis,
            fig_dnds_tru,
            go.Figure(), # do not display comparison plot

            # dndscv report figure style
            QQ_PLOT_RESULT_STYLE,
            DNDS_GLOBAL_RESULT_STYLE,
            DNDS_MIS_RESULT_STYLE,
            DNDS_TRUNC_RESULT_STYLE, 
            HIDE_PLOTS_STYLE,

            # dndscv report table style
            SHOW_TABLE_STYLE,
            HIDE_TABLE_STYLE,
            HIDE_TABLE_STYLE,
            HIDE_TABLE_STYLE,
            HIDE_TABLE_STYLE,
            
            dnd_radio_item_selection,
            RESULTS_DROPDOWN_OPTIONS,
            DNDS_RESULTS_DROPDOWM_LABEL,

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
    dropdown_menu_value,
    dnd_radio_item_selection           
):
    """ 
    """
    debugging= ""
    all_page_content = []
    figure1 = go.Figure()
    figure2 = go.Figure()
    figure3 = go.Figure()
    figure4 = go.Figure()
    comparison_fig = go.Figure()
    table1 = pd.DataFrame().to_dict('records')
    table2 = pd.DataFrame().to_dict('records')
    table3 = pd.DataFrame().to_dict('records')
    table4 = pd.DataFrame().to_dict('records')
    table5 = pd.DataFrame().to_dict('records')

    all_page_content = [
        table1,
        table2,
        table3,
        table4,
        table5,

        figure1,
        figure2, 
        figure3,
        figure4,
        comparison_fig,
        
        HIDE_PLOTS_STYLE, 
        HIDE_PLOTS_STYLE,
        HIDE_PLOTS_STYLE,
        HIDE_PLOTS_STYLE,
        DNDS_COMPARISON_STYLE,

        HIDE_TABLE_STYLE,
        SHOW_TABLE_STYLE,
        SHOW_TABLE_STYLE,
        SHOW_TABLE_STYLE,
        HIDE_TABLE_STYLE,

        dnd_radio_item_selection,
        COMPARISON_DROPDOWN_OPTIONS,
        DNDS_COMPARISON_DROPDOWM_LABEL,

        debugging
    ]

    return all_page_content

def gen_dndscv_summary_app_component(
    dnd_df_plot,
    dnd_df_merged,
    dnd_df_global,
    dropdown_menu_value,
    dnd_radio_item_selection       
):
    """
    """
    debugging= ""
    all_page_content = []
    figure1 = go.Figure()
    figure2 = go.Figure()
    figure3 = go.Figure()
    figure4 = go.Figure()
    figure5 = go.Figure()
    summary_table1 = dnd_df_plot.to_dict('records')
    table2 = pd.DataFrame().to_dict('records')
    table3 = pd.DataFrame().to_dict('records')
    table4 = pd.DataFrame().to_dict('records')
    table5 = dnd_df_merged.to_dict('records'),

    all_page_content = [
        # tables to display for the summary dNdScv report
        summary_table1,
        table2,
        table3,
        table4,
        table5,

        # no figures to display for the summary dNdScv report
        figure1,
        figure2, 
        figure3,
        figure4,
        figure5,

        HIDE_PLOTS_STYLE,
        HIDE_PLOTS_STYLE,
        HIDE_PLOTS_STYLE,
        HIDE_PLOTS_STYLE,
        HIDE_PLOTS_STYLE,

        SHOW_TABLE_STYLE,
        HIDE_TABLE_STYLE,
        HIDE_TABLE_STYLE,
        HIDE_TABLE_STYLE,
        SHOW_TABLE_STYLE,

        dnd_radio_item_selection,
        SUMMARY_DROPDOWN_OPTIONS,
        DNDS_SUMMARY_DROPDOWM_LABEL,

        debugging
    ]

    return all_page_content