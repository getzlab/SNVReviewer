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
from SNVReviewers.AppComponents.dNdScvComponentHelpers import gen_dndscv_results_app_component
from SNVReviewers.AppComponents.DIGAppComponent import DND_DATAFRAME_IDX, SNV_DATA_COLUMN_NAME

import pandas as pd
import numpy as np

DNDSCV_REPORT_VALUES = ["Results", "Comparison", "Summary"]

DNDSCV_REPORT_COLUMN_NAMES = ["RANK", "GENE", "N_SYN", "N_MIS", "N_NON", "N_SPL", "N_IND",
                              "dNdS_MIS", "dNdS_NON", "dNdS_SPL", "dNdS_IND", "PVAL_MIS",
                              "PVAL_TRUNC", "PVAL_IND", "PVAL", "FDR", "CGC", "PANCAN"]

def gen_dNdScv_app_component_data_callback(
    data: GenericData,
    idx,
    dnd_radio_item_selection,

):
    all_page_content = []
    dnd_df = data.df[SNV_DATA_COLUMN_NAME][0][DND_DATAFRAME_IDX]
    num_gene_values = 'All'
    figure1 = go.Figure()

    if dnd_radio_item_selection == "Results":
        all_page_content = gen_dndscv_results_app_component(
                dnd_df,
                num_gene_values,
                dnd_radio_item_selection
    )
    
    # elif dnd_radio_item_selection == "Comparison":
        # all_page_content = [
        #     figure1,
        #     dnd_radio_item_selection,
        #     COMPARISON_DROPDOWN_VALUES,
        # ]

    return all_page_content

def gen_dNdScv_app_component_layout():
    
    # table
    #
    return [        
            # Plotly Figure for the DIG Report
            html.Div([
                # radio button for selecting which type of report to display
                dbc.RadioItems(
                    options=[
                        {
                            "label": v, 
                            "value": v
                        } for v in DNDSCV_REPORT_VALUES
                    ],
                    value="Results",
                    id="dnds-report-type-radioitems",
                ),
                dbc.Row([
                    # insert the dropdown menus as columns inside this list for dbc.Row
                    dbc.Col([
                        # dropdown for selecting number of significant gene labels to display
                        dbc.Label("Select Number of Significant Gene Labels to Display:"),
                        dcc.Dropdown(
                        id='dnds-gene-dropdown',
                        options=[],
                        value='',
                        ),
                    ]),
                    
                ]),

                # REMOVE LATER!!!
                html.Div([
                    dbc.Label(id="dnds-special-text-output", children=""),
                ]),
                # Graphs above the coding region table
                dbc.Row([

                    
                    dbc.Col([
                        # creates the dig QQ plot
                        dcc.Graph(id='dnds-qq-graph', 
                                  figure={},
                                  style={"width":"1200px"} # increases the size of the plot
                                ), 
                    ]),
                    # # can hide the plots depending on the report
                    # dbc.Col([
                    #     # creates the dig QQ plot
                    #     dcc.Graph(id='dnd-volcano-graph', 
                    #               figure={},
                    #               style={"display":"none"} # hides the plot
                    #               ),
                    # ])
                ])
            ]),

            html.Div(
                [
                    # # displays the type of dig report you want displayed
                    # dbc.Row([
                    #     html.Div(
                    #         [
                    #             dbc.Label("dndSCV Report Table: "),
                    #             html.Label(children="Results", id="dnd-report-type-label"), # initialize label to empty string
                    #         ])
                    #     ]),                
                    
                # displays a table for the dig report
                html.Div(
                    children=[
                        html.H2('dNdScv Table'),
                        dash_table.DataTable(
                        id='dnds-report-table',
                        columns=[
                            {"name": i,
                                "id": i} for i in DNDSCV_REPORT_COLUMN_NAMES
                        ],
                        data=pd.DataFrame(columns=DNDSCV_REPORT_COLUMN_NAMES).to_dict(
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
                # Graphs below the coding region table
                dbc.Row([
                    dbc.Col([
                        # creates the dNdS ratio across all mutations plot
                        dcc.Graph(id='dnds-mutation-ratio-graph', 
                                  figure={},
                                  style={'display':'block', 'width':'600px'},
                                  ),
                    ]),
                    dbc.Col([
                        # creates the dNdS ratio of missense mutation plot
                        dcc.Graph(id='dnds-missense-graph', 
                                  figure={},
                                  style={'display':'block', 'width':'600px'},
                                  ),
                    ]),
                    dbc.Col([
                        # creates the dnds ratio of truncating mutations plot
                        dcc.Graph(id='dnds-truncating-graph', 
                                  figure={},
                                  style={'display':'block', 'width':'600px'}, 
                                  ),
                    ])
                ])
            ], 
        )
    ]

def gen_dnd_scv_app_component():
    
    return AppComponent(
        name='dNdScv Component',
        layout=gen_dNdScv_app_component_layout(),
        new_data_callback=gen_dNdScv_app_component_data_callback,
        internal_callback=gen_dNdScv_app_component_data_callback,
        callback_input=[
            Input('dnds-report-type-radioitems', 'value')
        ],
        callback_output=[
            Output('dnds-qq-graph', 'figure'),
            Output('dnds-special-text-output', 'children'),
            Output('dnds-gene-dropdown', 'options')
        ],
    )
