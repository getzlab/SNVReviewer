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

from SNVReviewers.AppComponents.DIGComponentHelpers import gen_combined_app_component, gen_coding_region_app_component, gen_non_coding_app_component
from SNVReviewers.AppComponents.utils import coding_region_mutation_type, combined_mutation_type, non_coding_region_burden_type

DIG_REPORT_COLUMN_NAMES = ["RANK", "GENE", "FDR", "PVAL", "PVAL_coding", "PVAL_promoter", "PVAL_5utr", 
                           "SIZE_coding", "SIZE_promoter", "SIZE_5utr", "SIZE_3utr", "CGC", "PANCAN"]
                    
DIG_REPORT_VALUES = ["Combined", "Coding region", "Promoter region", "5-prime UTRs", "3-prime UTRs", "Introns"]
# Different DIG Dataframes
DIG_COMBINED_DATAFRAME_IDX = 0
DIG_PRIME3_DATAFRAME_IDX = 1
DIG_PRIME5_DATAFRAME_IDX = 2
DIG_PROMOTER_DATAFRAME_IDX = 3
DIG_INTRON_DATAFRAME_IDX = 4

DND_PLOT_DATAFRAME_IDX = 5
DND_MERGED_DATAFRAME_IDX = 6
DND_GLOBAL_DATAFRAME_IDX = 7
MUTSIG_DATAFRAME_IDX = -1
SNV_DATA_COLUMN_NAME = "snv_data"
DIG_LABEL_IDX = 1

def gen_dig_app_component_data_internal_callback(
    data: GenericData,
    idx,
    # dig_label,
    dig_radio_item_selection,
    mutation_type,
    burden_type,
    p_val_type,
    display_toggle_value,
    display_label_value
):
    """
    
    """
    all_page_content = []
    dig_df = data.df[SNV_DATA_COLUMN_NAME][0][DIG_COMBINED_DATAFRAME_IDX].copy() # gets the dig report data for a specific cohort
    dig_prime3_df = data.df[SNV_DATA_COLUMN_NAME][0][DIG_PRIME3_DATAFRAME_IDX].copy()
    dig_prime5_df = data.df[SNV_DATA_COLUMN_NAME][0][DIG_PRIME5_DATAFRAME_IDX].copy()
    dig_promoter_df = data.df[SNV_DATA_COLUMN_NAME][0][DIG_PROMOTER_DATAFRAME_IDX].copy()
    dig_intron_df = data.df[SNV_DATA_COLUMN_NAME][0][DIG_INTRON_DATAFRAME_IDX].copy()

    if dig_radio_item_selection == 'Combined':
        
        # checking if you are changing to a new report type
        if mutation_type not in combined_mutation_type:
            # default values for the combined dig report 
            mutation_type = 'indels_snvs'
            burden_type = 'total'
            p_val_type = 'uniform_p_mid'

        all_page_content = gen_combined_app_component(
            dig_df,
            # dig_radio_item_selection, 
            mutation_type,
            burden_type,
            p_val_type,
            display_toggle_value,
            display_label_value
        )
        
    elif dig_radio_item_selection == 'Coding region':

        # checking if you are changing to a new report type
        if mutation_type not in coding_region_mutation_type:
            # default values for the coding region dig report
            mutation_type = 'indels_nonsynonymous_snvs'
            burden_type = 'total'
            p_val_type = 'uniform_p_mid'

        all_page_content = gen_coding_region_app_component(
            dig_df,
            # dig_radio_item_selection, 
            mutation_type,
            burden_type,
            p_val_type,
            display_toggle_value,
            display_label_value
        )

    # NEED TO FIGURE OUT HOW TO GET THE PROMOTER, 5 PRIME UTR REPORT TYPES WORKING
        # LOOK AT THE WORKFLOW THAT DAVID SENT
        # MIGHT NEED TO MAKE A SEPARATE LIST OF THE 5 DIFFERENT DATAFRAMES FOR ONE COHORT
            # MIGHT NEED TO HAVE DIFFERENT DATA PASSED INTO THE DIFFERENT gen_xxx_app_component functions
            # based on the data needed

    # THEN START WORKING ON THE mutsig component and dndscv component

    elif dig_radio_item_selection == "Promoter region":
        # checking if you are changing to a new report type
        if mutation_type not in combined_mutation_type:
            # default values for the 3 prime utr dig report
            mutation_type = 'indels_snvs'
            burden_type = 'total'
            p_val_type = 'uniform_p_mid'

        all_page_content = gen_non_coding_app_component(
            dig_promoter_df,
            # dig_radio_item_selection, 
            mutation_type,
            burden_type,
            p_val_type,
            display_toggle_value,
            display_label_value
        )

    elif dig_radio_item_selection == "5-prime UTRs":
        # checking if you are changing to a new report type
        if mutation_type not in combined_mutation_type:
            # default values for the 3 prime utr dig report
            mutation_type = 'indels_snvs'
            burden_type = 'total'
            p_val_type = 'uniform_p_mid'

        all_page_content = gen_non_coding_app_component(
            dig_prime5_df,
            # dig_radio_item_selection, 
            mutation_type,
            burden_type,
            p_val_type,
            display_toggle_value,
            display_label_value
        )

    elif dig_radio_item_selection == "3-prime UTRs":

        # checking if you are changing to a new report type
        if mutation_type not in combined_mutation_type:
            # default values for the 3 prime utr dig report
            mutation_type = 'indels_snvs'
            burden_type = 'total'
            p_val_type = 'uniform_p_mid'

        all_page_content = gen_non_coding_app_component(
            dig_prime3_df,
            # dig_radio_item_selection, 
            mutation_type,
            burden_type,
            p_val_type,
            display_toggle_value,
            display_label_value
        )

    elif dig_radio_item_selection == 'Introns':
        # checking if you are changing to a new report type
        if mutation_type not in combined_mutation_type:
            # default values for the 3 prime utr dig report
            mutation_type = 'indels_snvs'
            burden_type = 'total'
            p_val_type = 'uniform_p_mid'

        all_page_content = gen_non_coding_app_component(
            dig_intron_df,
            # dig_radio_item_selection, 
            mutation_type,
            burden_type,
            p_val_type,
            display_toggle_value,
            display_label_value
        )

    return all_page_content

def gen_dig_app_component_layout():
    """
    Generates the html layout for the DIG app component displaying the combined, coding region, 
    promoter, intron, 3 prime utr, and 5 prime utr DIG reports

    Parameters
    ==========
        None

    Return
    ======
        dash.html
            a plotly dash layout with a radio item for switching between the different report types, three dropdown menus corresponding
            to which mutation type/burden type/p value type you want to see displayed in the graphs and table, graphs for qq plot/volcano
            plot/fig mu plot/fig sigma plot/dnds fig plot, and warning text that displays if a specific combination of mutation and burden
            type are selected in the dropdown menu
    """
    return [
            html.Div(id='debugging-dig', children=''),
            # displays the interactive component to filter the samples displays based on their purity values
        
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
                            
                            # MIGHT BE GETTING RID OF THIS LATER!!!
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
                        dcc.Graph(id='dig-qq-graph', 
                                  figure={},
                                  style={"width":"1200px"} # increases the size of the plot
                                ), 
                    ]),
                    # can hide the plots depending on the report
                    dbc.Col([
                        # creates the dig QQ plot
                        dcc.Graph(id='dig-volcano-graph', 
                                  figure={},
                                  style={"display":"none"} # hides the plot
                                  ),
                    ])
                ])
            ]),

            html.Div(
                [
                    # # displays the type of dig report you want displayed
                    # dbc.Row([
                    #     html.Div(
                    #         [
                    #             dbc.Label("Dig Report Table: "),
                    #             html.Label(children="Combined", id="dig-report-type-label"), # initialize label to empty string
                    #         ])
                    #     ]),                
                    
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
                        dcc.Graph(id='dig-fig-mu-graph', 
                                  figure={},
                                  style={"display":"none"} # hides the plot
                                  ),
                    ]),
                    dbc.Col([
                        # creates the dig fig sigma plot
                        dcc.Graph(id='dig-fig-sigma-graph', 
                                  figure={},
                                  style={"display":"none"} # hides the plot
                                  ),
                    ]),
                    dbc.Col([
                        # creates the dig dnds fig plot
                        dcc.Graph(id='dig-dnds-fig-graph', 
                                  figure={},
                                  style={"display":"none"} # hides the plot
                                  ),
                    ])
                ])
            ], 
        )
    ]

def gen_dig_report_app_component():
    """
    Generates an AppComponent defining the interactive elements for viewing DIG reports
    
    Returns
    =======
    AnnoMate.AppComponent
        AppComponent defining the interactive elements for viewing DIG reports
    """
    
    return AppComponent(
        name='DIG Component',
        layout=gen_dig_app_component_layout(),
        new_data_callback=gen_dig_app_component_data_internal_callback,
        internal_callback=gen_dig_app_component_data_internal_callback,
        callback_input=[
            # Input('dig-report-type-label', 'children'),
            Input('dig-report-type-radioitems', 'value'), # mode value
            Input('dig-mutation-dropdown', 'value'),
            Input('dig-burden-dropdown', 'value'),
            Input('dig-p-value-dropdown', 'value'),
            Input('display-bounds-toggle-switch', 'on'),
            Input('display-labels-toggle-switch', 'on')
        ],

        callback_output=[
            Output('dig-report-table', 'data'),
            # Output('dig-report-type-label', 'children'),
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

            # changing the style of the plots to hide or display them and the size in which to display them
            Output('dig-qq-graph', 'style'),
            Output('dig-volcano-graph', 'style'),
            Output('dig-fig-mu-graph', 'style'),
            Output('dig-fig-sigma-graph', 'style'),
            Output('dig-dnds-fig-graph', 'style'),

            # warning message
            Output('special-text-output', 'children'),
            Output('debugging-dig', 'children')
        ],
    )
