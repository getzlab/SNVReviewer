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

def gen_dig_app_component_data_callback(
    data: GenericData,
    idx,
    dig_table_type_label
):
    # add a radio item that allows the user to be able to select which type of dig table they want to view
        # coding
        # 3-utr
        # 5-utr
        # non coding
        # combined
    
    

    return [
            data.df.to_dict('records'),
            "Coding"
            ]

def gen_dig_app_component_layout():
    
    # datatable
    #
    step_size = 5

    return [
            # displays the interactive component to filter the samples displays based on their purity values
            html.Div(
                [
                    # displays the type of dig report you want displayed
                    dbc.Row([
                        html.Div(
                            [
                                dbc.Label("Dig Report Table: "),
                                html.Label(children="", id="dig-report-type-label"), # initialize label to empty string
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
    
    return AppComponent(
        name='DIG Component',
        layout=gen_dig_app_component_layout(),
        new_data_callback=gen_dig_app_component_data_callback,
        internal_callback=gen_dig_app_component_data_callback,
        callback_input=[
            Input('dig-report-type-label', 'children')
            # Input('dig-report-coding-table', ''),
        ],
        callback_output=[
            Output('dig-report-coding-table', 'data'),
            Output('dig-report-type-label', 'children'),
        ],
    )
