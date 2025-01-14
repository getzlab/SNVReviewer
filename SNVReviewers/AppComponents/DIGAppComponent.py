from dash import dcc, html
from dash.dependencies import Input, Output
import dash_bootstrap_components as dbc
import plotly.graph_objects as go

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

from SNVReviewers.AppComponents.utils import generate_dig_report_plots, generate_plot_data

DIG_REPORT_COLUMN_NAMES = ["GENE", "CHROM", "GENE_LENGTH", "PVAL", "FDR", "OBS", "EXP", "MU", "SIGMA", ""]
                        # [
                        #     "GENE", "CHROM", "GENE_LENGTH", "R_SIZE", "R_OBS", "R_INDEL", 
                        #    "MU", "SIGMA", "ALPHA", "THETA", "MU_INDEL", "SIGMA_INDEL", "ALPHA_INDEL",
                        #    "THETA_INDEL", "FLAG", "Pi_SYN", "Pi_MIS", "Pi_NONS", "Pi_SPL", "Pi_TRUNC",
                        #    "Pi_NONSYN", "Pi_INDEL", "OBS_SYN", "OBS_MIS", "OBS_NONS", "OBS_SPL", "OBS_INDEL",
                        #    "OBS_TRUNC", "OBS_NONSYN", "N_SAMP_SYN",	"N_SAMP_MIS", "N_SAMP_NONS", "N_SAMP_SPL",
                        #    "N_SAMP_TRUNC", "N_SAMP_NONSYN",	"N_SAMP_INDEL",	"EXP_SYN", "EXP_MIS", "EXP_NONS",	
                        #    "EXP_SPL", "EXP_TRUNC", "EXP_NONSYN", "PVAL_SYN_BURDEN",	"PVAL_MIS_BURDEN",	"PVAL_NONS_BURDEN",	
                        #    "PVAL_SPL_BURDEN", "PVAL_TRUNC_BURDEN",	"PVAL_NONSYN_BURDEN", "PVAL_SYN_BURDEN_SAMPLE",
                        #    "PVAL_MIS_BURDEN_SAMPLE", "PVAL_NONS_BURDEN_SAMPLE",	"PVAL_SPL_BURDEN_SAMPLE", "PVAL_TRUNC_BURDEN_SAMPLE",	
                        #    "PVAL_NONSYN_BURDEN_SAMPLE",	"EXP_INDEL", "PVAL_INDEL_BURDEN", "PVAL_MUT_BURDEN"
                        # ]
DIG_REPORT_VALUES = ["Combined", "Coding regions", "Promoter regions", "5-prime UTRs", "3-prime UTRs"]
DIG_DATAFRAME_IDX = 0
DND_DATAFRAME_IDX = 1
MUTSIG_DATAFRAME_IDX = 2
SNV_DATA_COLUMN_NAME = "snv_data"
DIG_LABEL_IDX = 1
# MAKE RADIO ITEM FOR THE COMBINED REPORT
    # RADIO ITEM WILL SWITCH TO THE TXT FILE THAT IS USED TO CREATE PLOTS/TABLE
# MAKE RADIO ITEM FOR THE CODING REGION REPORT

def gen_dig_app_component_data_internal_callback(
    data: GenericData,
    idx,
    dig_label,
    dig_type_selection
):
    """
    
    """
    # add a radio item that allows the user to be able to select which type of dig table they want to view
        # coding
        # 3-utr
        # 5-utr
        # non coding
        # combined
    # print(type)
    dig_df = data.df[SNV_DATA_COLUMN_NAME][0][DIG_DATAFRAME_IDX] # gets the dig report data for a specific cohort
    dnd_df = data.df[SNV_DATA_COLUMN_NAME][0][DND_DATAFRAME_IDX]
    mutsig_df = data.df[SNV_DATA_COLUMN_NAME][0][MUTSIG_DATAFRAME_IDX]
    debugging_component = "" + str(type(data.df["snv_data"]))
    dig_output_type = "Combined"
    display_bounds = True # whether to display the bounds
    
    dir_output = "./example_notebooks/data/dig_data"
    dig_figure = go.Figure()
    qq_fig = go.Figure()

    # do not need to run this because the dig report data is already generated
    # generate_dig_report(
    #     dig_df,
    #     dig_output_type,
    #     dir_output,
    # )
    qq_fig, table_fig = generate_dig_report_plots(dig_df)

    # ONLY GETTING THE FIRST 500 ROWS OF DATA TO DISPLAY IN THE TABLE   
    dig_df = dig_df[:500]


    # wrap up the precalled purity 

    return [
            dig_df.to_dict('records'),
            dig_type_selection,
            qq_fig,
            debugging_component,
            ]

def gen_dig_app_component_data_external_callback(
    data: GenericData,
    idx,
    dig_label,
    dig_type_selection
):
    """
    """
    # add a radio item that allows the user to be able to select which type of dig table they want to view
        # coding
        # 3-utr
        # 5-utr
        # non coding
        # combined
    # print(type)
    output = gen_dig_app_component_data_internal_callback(
                data,
                idx,
                dig_label,
                dig_type_selection
            )
    
    # output[DIG_LABEL_IDX] = dig_label

    return output

def gen_dig_app_component_layout():
    """
    
    """
    
    # datatable
    #
    step_size = 5

    return [
            # displays the interactive component to filter the samples displays based on their purity values
            html.Div([
                dbc.Label(children="Debugging Stuff!!!", id="debugging"),    
            ]),

            # Plotly Figure for the DIG Report
            html.Div([

                # radio button for selecting which type of report to display
                dbc.RadioItems(
                    options=[
                        {
                            "label": v, 
                            "value": v
                        } for v in DIG_REPORT_VALUES
                    ],
                    value="Combined",
                    id="dig-report-type-radioitems",
                ),

                # dropdown for selecting a mutation type
                dbc.Label("Select Mutation Type"),
                dcc.Dropdown(
                id='dig-mutation-dropdown',
                options=[
                    {'label': 'Indels + SNVs', 'value': 'indels_snvs'},
                    {'label': 'Indels', 'value': 'indels'},
                    {'label': 'SNVs', 'value': 'snvs'}
                ]),


                # dropdown for selecting burden type
                dbc.Label("Select Burden Type"),
                dcc.Dropdown(
                id='dig-burden-dropdown',
                options=[
                    {'label': 'Total', 'value': 'total'},
                    {'label': 'Sample-wise', 'value': 'sample_wise'}
                ]),

                # dropdown for selecting burden type
                dbc.Label("P-value Type"),
                dcc.Dropdown(
                id='dig-p-value-dropdown',
                options=[
                    {'label': 'Uniform P-mid', 'value': 'uniform_p_mid'},
                    {'label': 'P-mid', 'value': 'p_mid'},
                ]),

                # creates the dig QQ plot
                dcc.Graph(id='dig-qq-graph', figure={}),
            ]),

            html.Div(
                [
                    # displays the type of dig report you want displayed
                    dbc.Row([
                        html.Div(
                            [
                                dbc.Label("Dig Report Table: "),
                                html.Label(children="Combined", id="dig-report-type-label"), # initialize label to empty string
                            ])
                        ]),                
                    
                # displays a table for the dig report
                html.Div(
                    children=[
                        html.H2('DIG Table'),
                        dash_table.DataTable(
                        id='dig-report-coding-table',
                        columns=[
                            {"name": i,
                                "id": i} for i in DIG_REPORT_COLUMN_NAMES
                        ],
                        data=pd.DataFrame(columns=DIG_REPORT_COLUMN_NAMES).to_dict(
                            'records'),
                        editable=False,
                        filter_action="native",
                        sort_action="native",
                        sort_mode="multi",
                        row_selectable="single",
                        row_deletable=False,
                        selected_columns=[],
                        selected_rows=[0],
                        page_action="native",
                        page_current=0,
                        page_size=5),
                    ]
                )
            ], 
        )
    ]

def gen_dig_report_app_component():
    """
    
    """
    
    return AppComponent(
        name='DIG Component',
        layout=gen_dig_app_component_layout(),
        new_data_callback=gen_dig_app_component_data_internal_callback,
        internal_callback=gen_dig_app_component_data_external_callback,
        callback_input=[
            Input('dig-report-type-label', 'children'),
            Input('dig-report-type-radioitems', 'value'),
        ],

        callback_output=[
            Output('dig-report-coding-table', 'data'),
            Output('dig-report-type-label', 'children'),
            Output('dig-qq-graph', 'figure'),
            Output('debugging', 'children'),
        ],
    )
