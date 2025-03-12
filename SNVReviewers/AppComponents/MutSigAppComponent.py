from dash import dcc, html
import dash_bootstrap_components as dbc
from AnnoMate.ReviewDataApp import AppComponent
from AnnoMate.DataTypes.GenericData import GenericData
from dash import dcc, html, dash_table
from dash.dependencies import Input, Output
import dash_bootstrap_components as dbc

from AnnoMate.ReviewDataApp import AppComponent
from AnnoMate.DataTypes.GenericData import GenericData
from SNVReviewers.AppComponents.DIGAppComponent import MUTSIG_DATAFRAME_IDX
from SNVReviewers.AppComponents.DIGAppComponent import SNV_DATA_COLUMN_NAME
from SNVReviewers.AppComponents.dNdScvComponentHelpers import RESULTS_DROPDOWN_OPTIONS
from SNVReviewers.AppComponents.MutSigComponentHelpers import gen_mutsig_results_app_component, MUTSIG_REPORT_COLUMN_NAMES

import pandas as pd
import numpy as np



def gen_mutsig_app_component_data_callback(
    data: GenericData,
    idx,
    mutsig_dropdown_value
):
    """
    """
    all_page_content = []
    mutsig_df = data.df[SNV_DATA_COLUMN_NAME][0][MUTSIG_DATAFRAME_IDX]
    mustsig_dropwdown_values = [option["value"] for option in RESULTS_DROPDOWN_OPTIONS]

    if mutsig_dropdown_value not in mustsig_dropwdown_values:
        mutsig_dropdown_value = "all"

    all_page_content = gen_mutsig_results_app_component(mutsig_df, mutsig_dropdown_value)
    
    return all_page_content

def gen_mutsig_app_component_layout():
    """ 
    """
    
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
                            dbc.Label(id="mutsig-dropdowm-label", children="Number of Significant Gene Labels to Display: "),
                        ]),
                        dcc.Dropdown(
                        id='mutsig-gene-dropdown',
                        options=RESULTS_DROPDOWN_OPTIONS,
                        value='',
                        ),
                    ]),
                    
                ]),

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
                            page_size=15,
                        
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
    """
    """
    
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
            Output('mutsig-qq-graph', 'figure'),
            Output('mutsig-gene-dropdown', 'value'),
            
            # REMOVE LATER!!
            Output('mutsig-debugging', 'children'),
        ],
    )
