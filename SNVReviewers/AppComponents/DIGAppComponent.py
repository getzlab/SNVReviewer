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

from SNVReviewers.AppComponents.DIGAppComponent_Combined import gen_dig_combined_app_component_data_internal_callback, gen_combined_dig_report_app_component, COMBINED_MUT_DROPDOWN, COMBINED_BUR_DROPDOWN, COMBINED_SCATTER_DROPDOWN
from SNVReviewers.AppComponents.DIGAppComponent_CodingRegion import gen_dig_coding_region_app_component_data_internal_callback, gen_dig_coding_region_app_component, gen_dig_app_coding_region_component_layout, DIG_CODING_REGION_REPORT_COLUMN_NAMES
from SNVReviewers.AppComponents.utils import generate_combined_dig_report_plots, generate_coding_region_report, coding_region_mutation_type, coding_region_burden_type, combined_mutation_type, scatterpoint_type

DIG_REPORT_COLUMN_NAMES = ["RANK", "GENE", "FDR", "PVAL", "PVAL_coding", "PVAL_promoter", "PVAL_5utr", 
                           "SIZE_coding", "SIZE_promoter", "SIZE_5utr", "SIZE_3utr", "CGC", "PANCAN"]
                    
DIG_REPORT_VALUES = ["Combined", "Coding region", "Promoter region", "5-prime UTRs", "3-prime UTRs"]
DIG_DATAFRAME_IDX = 0
DND_DATAFRAME_IDX = 1
MUTSIG_DATAFRAME_IDX = 2
SNV_DATA_COLUMN_NAME = "snv_data"
DIG_LABEL_IDX = 1

def gen_dig_app_component_data_internal_callback(
    data: GenericData,
    idx,
    dig_label,
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

    if dig_type_selection == 'Combined':
        
        # checking if you are changing to a new report type
        if mutation_type not in combined_mutation_type:
            mutation_type = 'indels_snvs'
            burden_type = 'total'
            p_val_type = 'uniform_p_mid'


        # PUT THESE FUNCTIONS INTO A HELPER FILE AND IMPORT FOR A
        all_page_content = gen_dig_combined_app_component_data_internal_callback(
            data,
            idx,
            dig_label,
            dig_type_selection, 
            mutation_type,
            burden_type,
            p_val_type,
            display_toggle_value,
            display_label_value
        )
        
    elif dig_type_selection == 'Coding region':
        # checking if you are changing to a new report type
        if mutation_type not in coding_region_mutation_type:
            mutation_type = 'indels_nonsynonymous_snvs'
            burden_type = 'total'
            p_val_type = 'uniform_p_mid'

        all_page_content = gen_dig_coding_region_app_component_data_internal_callback(
            data,
            idx,
            dig_label,
            dig_type_selection, 
            mutation_type,
            burden_type,
            p_val_type,
            display_toggle_value,
            display_label_value
        )

    elif dig_type_selection == "Promoter region":
        raise NotImplementedError
        print("I am in the promoter region")

    elif dig_type_selection == "5-prime UTRs":
        raise NotImplementedError
        print("I am in the 5 prime utrs region")

    elif dig_type_selection == "3-prime UTRs":
        raise NotImplementedError
        print("I am in the 3 prime utrs region")  

    return all_page_content

def gen_dig_app_component_data_external_callback(
    data: GenericData,
    idx,
    dig_label,
    dig_type_selection,
    mutation_type,
    burden_type,
    p_val_type,
    display_toggle_value, 
    display_label_value
):
    """

    """
    output = gen_dig_app_component_data_internal_callback(
                data,
                idx,
                dig_label,
                dig_type_selection, 
                mutation_type,
                burden_type,
                p_val_type,
                display_toggle_value,
                display_label_value
            )
    
    return output



def gen_dig_app_component_layout():
    """
    
    """
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
                dbc.Row([
                    # insert the dropdown menus as columns inside this list for dbc.Row
                    dbc.Col([
                        # dropdown for selecting a mutation type
                        dbc.Label("Select Mutation Type"),
                        dcc.Dropdown(
                        id='dig-mutation-dropdown',
                        options=[],
                        value='',
                        ),
                    ]),
                    dbc.Col([
                        # dropdown for selecting burden type
                        dbc.Label("Select Burden Type"),
                        dcc.Dropdown(
                        id='dig-burden-dropdown',
                        options=[],
                        value=''
                        ),
                    ]),
                    dbc.Col([
                        # dropdown for selecting burden type
                        dbc.Label("P-value Type"),
                        dcc.Dropdown(
                        id='dig-p-value-dropdown',
                        options=[],
                        value=''),
                    ])
                ]),
                html.Div([
                    dbc.Row([
                        dbc.Col([
                            # makes a toggle component
                            daq.BooleanSwitch(
                            id='display-bounds-toggle-switch',
                            label='Display Bounds',
                            on=False),
                        ]),
                        dbc.Col([
                            # makes a toggle component
                            daq.BooleanSwitch(
                            id='display-labels-toggle-switch',
                            label='Display Labels',
                            on=False),
                        ])
                    ])
                ]),

                html.Div([
                    dbc.Label(id="special-text-output", children=""),
                ]),
                # Graphs above the coding region table
                dbc.Row([
                    dbc.Col([
                        # creates the dig QQ plot
                        dcc.Graph(id='dig-qq-graph', figure={}),
                    ]),
                    dbc.Col([
                        # creates the dig QQ plot
                        dcc.Graph(id='dig-volcano-graph', figure={}),
                    ])
                ])
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
                        id='dig-report-table',
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
                        # dcc.Graph(id='dig-fig-mu-coding-region-graph', figure={}),
                        dcc.Graph(id='dig-fig-mu-graph', figure={}),
                    ]),
                    dbc.Col([
                        # creates the dig fig sigma plot
                        # dcc.Graph(id='dig-fig-sigma-coding-region-graph', figure={}),
                        dcc.Graph(id='dig-fig-sigma-graph', figure={}),
                    ]),
                    dbc.Col([
                        # creates the dig dnds fig plot
                        # dcc.Graph(id='dig-dnds-fig-coding-region-graph', figure={}),
                        dcc.Graph(id='dig-dnds-fig-graph', figure={}),
                    ])
                ])
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
        internal_callback=gen_dig_app_component_data_internal_callback,

        # MAYBE USE THE DROPDOWN MENUS AS INPUT AND OUTPUT (For changing which columns get accessed in the table!!)
            # SO THAT WHEN THE DROP DOWN CHANGES 

        # Need to do a conditional statement for whether to use these values for input or other report types for input
        callback_input=[
            Input('dig-report-type-label', 'children'),
            Input('dig-report-type-radioitems', 'value'), # mode value
            Input('dig-mutation-dropdown', 'value'),
            Input('dig-burden-dropdown', 'value'),
            Input('dig-p-value-dropdown', 'value'),
            Input('display-bounds-toggle-switch', 'on'),
            Input('display-labels-toggle-switch', 'on')
        ],

        callback_output=[
            # Output('final_container', 'children')
            # INSERT THE REMAINING POSSIBLE GRAPHS REQUIRED FOR MAKING ALL THE OTHER REPORT TYPES FOR THE DIG APP COMPONENT 
            # Output('dig-report-coding-table', 'data'),
            # Output('dig-report-type-label', 'children'),
            # Output('dig-qq-graph', 'figure'),
            # Output('dig-report-coding-table', 'columns'),
            # Output('dig-mutation-dropdown', 'value'),
            # Output('dig-burden-dropdown', 'value'),
            # Output('dig-p-value-dropdown', 'value'),
            # Output('special-text-output', 'children'),
            # Output('debugging', 'children'),

# trigger id, look into this
            Output('dig-report-table', 'data'),
            Output('dig-report-type-label', 'children'),
            Output('dig-volcano-graph', 'figure'),
            Output('dig-qq-graph', 'figure'),
            Output('dig-fig-mu-graph', 'figure'),
            Output('dig-fig-sigma-graph', 'figure'),
            Output('dig-dnds-fig-graph', 'figure'),
            Output('dig-report-table', 'columns'),
            Output('dig-mutation-dropdown', 'value'),
            Output('dig-burden-dropdown', 'value'),
            Output('dig-p-value-dropdown', 'value'),

            # returns the dropdown options
            Output('dig-mutation-dropdown', 'options'),
            Output('dig-burden-dropdown', 'options'),
            Output('dig-p-value-dropdown', 'options'),

            Output('special-text-output', 'children'),
            Output('debugging', 'children'),
        ],
    )
