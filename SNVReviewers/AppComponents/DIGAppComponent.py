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

from SNVReviewers.AppComponents.utils import generate_dig_report_plots, generate_plot_data, mut_type, scatterpoint_type

DIG_REPORT_COLUMN_NAMES = ["RANK", "GENE", "FDR", "PVAL", "PVAL_coding", "PVAL_promoter", "PVAL_5utr", 
                           # NEED TO ADD ON recalc, unif if the dropdown menu value is uniform or p-mid
                           # PVAL_coding_SNV -> PVAL_coding_SNV_recalc
                           "SIZE_coding", "SIZE_promoter", "SIZE_5utr", "SIZE_3utr", "CGC", "PANCAN"]
                    
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
    dig_type_selection,
    mutation_type,
    burden_type,
    p_val_type,
    display_toggle_value
):
    """
    
    """

    if mutation_type == "":
        mutation_type = 'indels_snvs'
        burden_type = 'total'
        p_val_type = 'uniform_p_mid'

    dig_df = data.df[SNV_DATA_COLUMN_NAME][0][DIG_DATAFRAME_IDX].copy() # gets the dig report data for a specific cohort
    dnd_df = data.df[SNV_DATA_COLUMN_NAME][0][DND_DATAFRAME_IDX].copy()
    mutsig_df = data.df[SNV_DATA_COLUMN_NAME][0][MUTSIG_DATAFRAME_IDX].copy()

    debugging_component = "" #+ str(type(data.df["snv_data"]))
    dig_output_type = "Combined"
    display_bounds = False # whether to display the bounds
    
    dir_output = "./example_notebooks/data/dig_data"
    dig_figure = go.Figure()
    qq_fig = go.Figure()

    # MAKE SURE TO ADD PLOTLY DASH COMPONENT THAT WILL UPDATE THIS VALUE
    display_bounds = 'No'
    display_labels_key = 'No'

    # ADD IN A SELECTION TOOL FOR DISPLAYING BOUNDS, KEYS FOR SEEING BOUNDS IS 'Yes' and 'No'
    qq_fig, table_fig, text_special = generate_dig_report_plots(dig_df, mutation_type, burden_type, display_bounds, display_labels_key, p_val_type)
    dig_data_columns = []
    
    for clm_nm in DIG_REPORT_COLUMN_NAMES:
        dig_data_clm_dict = {}
        dig_data_clm_dict["name"] = clm_nm
        dig_data_clm_dict["id"] = clm_nm

        if clm_nm == 'FDR' or "PVAL" in clm_nm:
            # gets the column name corresponding to the mutation type and the scatterpoint type
            dig_data_clm_dict["id"] = clm_nm + "_"+ mut_type[mutation_type] + "_" + scatterpoint_type[p_val_type]

        # if display_bounds:
        #     dig_data_clm_dict["id"] += "_lower"
        #     dig_data_columns.append(di g_data_clm_dict) # appends the lower bound column info
        #     dig_data_clm_dict = {}

        # debugging_component += "name: " + dig_data_clm_dict["name"] + " id: " + dig_data_clm_dict["id"] + "; "
        dig_data_columns.append(dig_data_clm_dict)
    
    for column_dict in dig_data_columns:
        # gets the dig data column name
        column = column_dict["id"]
        # debugging_component += "column name: " + column +"; "

        # skip the rank column
        if column == 'RANK':
            continue

        # rounds all the values in the FDR and PVAL columns to 3 decimal places
        if 'FDR' in column or 'PVAL' in column:
            dig_df[column] = np.round(dig_df[column], decimals=3)
        
    # get the coding region working plots working!!
    # get the display bounds selection tool working 

    # ONLY GETTING THE FIRST 500 ROWS OF DATA TO DISPLAY IN THE TABLE
    # REMOVE THE DEBUGGING LATER!!!   
    dig_df = dig_df[:500]

    # wrap up the precalled purity 

    return [
            dig_df.to_dict('records'),
            dig_type_selection,
            qq_fig,
            dig_data_columns,
            mutation_type,
            burden_type,
            p_val_type,
            text_special,
            debugging_component,
        ]

def gen_dig_app_component_data_external_callback(
    data: GenericData,
    idx,
    dig_label,
    dig_type_selection,
    mutation_type,
    burden_type,
    p_val_type,
    display_toggle_value
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
                display_toggle_value
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

                # dbc.Row([
                #     dbc.Col([
                #         dbc.RadioItems(
                #         options=[
                #             {
                #                 "label": val, 
                #                 "value": val
                #             }
                #         ],
                #         value="Combined",
                #         id=f"dig-report-type-radioitems-{DIG_REPORT_VALUES[idx]}",
                #         ),
                #     ]) for idx, val in enumerate(DIG_REPORT_VALUES)
                # ]),
                
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

                dbc.Row([# insert the dropdown menus as columns inside this list for dbc.Row
                    dbc.Col([
                        # dropdown for selecting a mutation type
                        dbc.Label("Select Mutation Type"),
                        dcc.Dropdown(
                        id='dig-mutation-dropdown',
                        options=[
                            {'label': 'Indels + SNVs', 'value': 'indels_snvs'},
                            {'label': 'Indels', 'value': 'indels'},
                            {'label': 'SNVs', 'value': 'snvs'}
                        ],
                        value='',
                        ),
                    ]),
                    dbc.Col([
                        # dropdown for selecting burden type
                        dbc.Label("Select Burden Type"),
                        dcc.Dropdown(
                        id='dig-burden-dropdown',
                        options=[
                            {'label': 'Total', 'value': 'total'},
                            {'label': 'Sample-wise', 'value': 'sample_wise'}
                        ],
                        value=''
                        ),
                    ]),
                    dbc.Col([
                        # dropdown for selecting burden type
                        dbc.Label("P-value Type"),
                        dcc.Dropdown(
                        id='dig-p-value-dropdown',
                        options=[
                            {'label': 'Uniform P-mid', 'value': 'uniform_p_mid'},
                            {'label': 'P-mid', 'value': 'p_mid'},
                        ],
                        value=''),
                    ])
                ]),
                html.Div([
                    dbc.Row([
                        dbc.Col([
                            html.Label(children="Display Bounds"), 
                        ]),
                        dbc.Col([
                            # makes a toggle component
                            daq.ToggleSwitch(
                            id='display-bounds-toggle-switch',
                            value=False),
                        ])
                    ])
                ]),

                html.Div([
                    dbc.Label(id="special-text-output", children=""),
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
                        page_size=5,
                        # changing the width of the data table to 
                        style_table={
                            'width': '100%',  # Make the table width responsive
                            'maxWidth': '100%',  # Ensure it doesn’t go beyond the screen width
                            'overflowX': 'auto',  # Allow horizontal scroll if necessary
                        },
                       
                        ),
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
            Input('dig-mutation-dropdown', 'value'),
            Input('dig-burden-dropdown', 'value'),
            Input('dig-p-value-dropdown', 'value'),
            Input('display-bounds-toggle-switch', 'value')
        ],

        callback_output=[
            Output('dig-report-coding-table', 'data'),
            Output('dig-report-type-label', 'children'),
            Output('dig-qq-graph', 'figure'),
            Output('dig-report-coding-table', 'columns'),
            Output('dig-mutation-dropdown', 'value'),
            Output('dig-burden-dropdown', 'value'),
            Output('dig-p-value-dropdown', 'value'),
            Output('special-text-output', 'children'),
            Output('debugging', 'children'),
        ],
    )
