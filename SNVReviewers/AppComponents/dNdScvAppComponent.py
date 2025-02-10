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

import pandas as pd
import numpy as np

DNDSCV_REPORT_VALUES = ["Results", "Comparison", "Summary"]
NUM_GENE_LABELS = ["All", "First 30", "First 20", "First 10", "None"]
DNDSCV_REPORT_COLUMN_NAMES = ["GENE", "N_SYN", "N_MIS"] # finish the rest of the column name

def gen_dNdScv_app_component_data_callback(
    data: GenericData,
    idx,
    dnd_radio_item_selection,

):
    all_page_content = []

    figure1 = go.Figure()
    all_page_content = [
        figure1,
        dnd_radio_item_selection,
        NUM_GENE_LABELS
    ]

    return all_page_content

def gen_dNdScv_app_component_layout():
    
    # table
    #
    return [
            # displays the interactive component to filter the samples displays based on their purity values
        
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
                    id="dnd-report-type-radioitems",
                ),
                dbc.Row([
                    # insert the dropdown menus as columns inside this list for dbc.Row
                    dbc.Col([
                        # dropdown for selecting number of significant gene labels to display
                        dbc.Label("Select Number of Significant Gene Labels to Display:"),
                        dcc.Dropdown(
                        id='dnd-gene-dropdown',
                        options=[],
                        value='',
                        ),
                    ]),
                    
                ]),
                html.Div([
                    dbc.Label(id="dnd-special-text-output", children=""),
                ]),
                # Graphs above the coding region table
                dbc.Row([

                    
                    dbc.Col([
                        # creates the dig QQ plot
                        dcc.Graph(id='dnd-qq-graph', 
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
                        html.H2('DND SCV Table'),
                        dash_table.DataTable(
                        id='dnd-report-table',
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
                        # creates the dig fig mu plot
                        dcc.Graph(id='dnd-mutation-ratio-graph', 
                                  figure={},
                                #   style={"display":"none"} # hides the plot
                                  ),
                    ]),
                    dbc.Col([
                        # creates the dig fig sigma plot
                        dcc.Graph(id='dnd-missense-graph', 
                                  figure={},
                                #   style={"display":"none"} # hides the plot
                                  ),
                    ]),
                    dbc.Col([
                        # creates the dig dnds fig plot
                        dcc.Graph(id='dnd-truncating-graph', 
                                  figure={},
                                #   style={"display":"none"} # hides the plot
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
            Input('dnd-report-type-radioitems', 'value')
        ],
        callback_output=[
            Output('dnd-qq-graph', 'figure'),
            Output('dnd-special-text-output', 'children'),
            Output('dnd-gene-dropdown', 'options')
        ],
    )
