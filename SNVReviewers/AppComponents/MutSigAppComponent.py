from dash import dcc, html
import dash_bootstrap_components as dbc
import plotly.graph_objects as go
import dash_daq as daq

from AnnoMate.ReviewDataApp import AppComponent
from AnnoMate.DataTypes.GenericData import GenericData

from dash import dcc, html, dash_table
from dash.dependencies import Input, Output, State
import dash_bootstrap_components as dbc
import plotly.graph_objects as go

from AnnoMate.Data import Data, DataAnnotation
from AnnoMate.ReviewDataApp import ReviewDataApp, AppComponent
from AnnoMate.DataTypes.GenericData import GenericData
from SNVReviewers.AppComponents.dNdScvComponentHelpers import gen_dndscv_results_app_component, gen_dndscv_comparison_app_component, gen_dndscv_summary_app_component
from SNVReviewers.AppComponents.DIGAppComponent import DND_PLOT_DATAFRAME_IDX, DND_MERGED_DATAFRAME_IDX, DND_GLOBAL_DATAFRAME_IDX, DND_COMPARISON_DATAFRAME_IDX
from SNVReviewers.AppComponents.DIGAppComponent import SNV_DATA_COLUMN_NAME
from SNVReviewers.AppComponents.dNdScvComponentHelpers import RESULTS_DROPDOWN_OPTIONS, COMPARISON_DROPDOWN_OPTIONS, SUMMARY_DROPDOWN_OPTIONS

import pandas as pd
import numpy as np

MUTSIG_REPORT_COLUMN_NAMES = ["RANK", "GENE", "NNEI", "NNCD", "NSIL", "NMIS", "NSTP", "NSPL", "NIND"
                              #, # FINISH PUTTING THE REST OF THE COLUMN NAMES LATER!!
                              ]

def gen_mutsig_app_component_data_callback(
    data: GenericData,
    idx,
    mutsig_dropdown_value
):
    data = pd.DataFrame().to_dict("records")
    return[ 
        data
    ]

def gen_mutsig_app_component_layout():
    
    return [   
            # REMOVE LATER!!!
            html.Div([
                dbc.Label(id="mutsig-debugging", children=""),
            ]),

            # Plotly Figure for the DIG Report
            html.Div([
                dbc.Row([
                    # insert the dropdown menus as columns inside this list for dbc.Row
                    dbc.Col([
                        # dropdown for selecting number of significant gene labels to display
                        html.Div([
                            dbc.Label(id="mutsig-dropdowm-label", children=""),
                        ]),
                        # dbc.Label('dnds-dropdowm-label', children=""),
                        dcc.Dropdown(
                        id='mutsig-gene-dropdown',
                        options=[],
                        value='',
                        ),
                    ]),
                    
                ]),

                # REMOVE LATER!!!
                html.Div([
                    dbc.Label(id="mutsig-special-text-output", children=""),
                ]),
                # Graphs above the coding region table
                dbc.Row([
                    dbc.Col([
                        # creates the mutsig QQ plot
                        dcc.Graph(id='mutsig-qq-graph', 
                                  figure={},
                                  style={"width":"1200px"} # increases the size of the plot
                                ), 
                    ]),
                ])
            ]),

            html.Div(
                [
                    # displays a table for the dig report
                    html.Div(
                        children=[
                            # html.H2('dNdScv Table'),
                            dash_table.DataTable(
                            id='mutsig-report-table',
                            columns=[
                                {"name": i,
                                    "id": i} for i in MUTSIG_REPORT_COLUMN_NAMES
                            ],
                            data=pd.DataFrame(columns=MUTSIG_REPORT_COLUMN_NAMES).to_dict(
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
                            page_size=5,
                        
                            # changing the width of the data table to 
                            style_table={
                                'width': '100%',  # Make the table width responsive
                                'maxWidth': '100%',  # Ensure it doesn’t go beyond the screen width
                                'overflowX': 'auto',  # Allow horizontal scroll if necessary
                                },
                            ),
                        ]
                    ),
                ], 
            )
        ]

def gen_mutsig_app_component():
    
    return AppComponent(
        name='MutSig Component',
        layout=gen_mutsig_app_component_layout(),
        new_data_callback=gen_mutsig_app_component_data_callback,
        internal_callback=gen_mutsig_app_component_data_callback,
        callback_input=[
            Input('mutsig-gene-dropdown', 'value')
        ],
        callback_output=[
            Output('mutsig-report-table', 'data'),
        ],
    )
