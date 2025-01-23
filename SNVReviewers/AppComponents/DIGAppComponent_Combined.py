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

from SNVReviewers.AppComponents.utils import generate_combined_dig_report_plots, generate_coding_region_report, coding_region_mutation_type, coding_region_burden_type, combined_mutation_type, combined_burden_type, scatterpoint_type

DIG_REPORT_COLUMN_NAMES = ["RANK", "GENE", "FDR", "PVAL", "PVAL_coding", "PVAL_promoter", "PVAL_5utr", 
                           "SIZE_coding", "SIZE_promoter", "SIZE_5utr", "SIZE_3utr", "CGC", "PANCAN"
                            ]
CODING_REGION_DIG_REPORT_COLUMN_NAMES = ["RANK", "GENE", "CHROM", "LENGTH", "FDR", "PVAL", "OBS", "EXP",
                                         "MU", "SIGMA", "dNdS_OBS", "dNdS_EXP", "FLAG", "CGC", "PANCAN"
                                        ]
                    
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
    display_toggle_value,
    display_label_value
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

    # if dig_label == "Combined":
    if combined_burden_type[burden_type] == "":
        dig_df = dig_df.sort_values(by='PVAL' + "_" + combined_mutation_type[mutation_type] + "_" + scatterpoint_type[p_val_type])
    else:
        dig_df = dig_df.sort_values(by='PVAL'+ "_"+ combined_mutation_type[mutation_type] + "_" + combined_burden_type[burden_type] + "_" + scatterpoint_type[p_val_type])
    dig_df['RANK'] = np.array([i+1 for i in range(len(dig_df))])

    debugging_component = ""
    dig_output_type = "Combined"
    
    dir_output = "./example_notebooks/data/dig_data"
    dig_figure = go.Figure()
    qq_fig = go.Figure()

    # checks if the user wants to display the bounds on the dig report plot
    if display_toggle_value:
        display_bounds = 'Yes'

    # defaults to not displaying lower/upper bounds on the dig report plot
    else:
        display_bounds = 'No'
    
    if display_label_value:
        display_labels_key = 'Yes'
    
    else:
        display_labels_key = 'No'

    qq_fig, table_fig, text_special = generate_combined_dig_report_plots(dig_df, mutation_type, 
                                                                         burden_type, display_bounds, 
                                                                         display_labels_key, p_val_type)
    

    # MODIFY THE LAYOUT HERE, CALL THE CODING/COMBINED APP LAYOUT FUNCTION
    if dig_label == "Combined":
        children_layout = gen_dig_app_combined_component_layout()

    # # MAKE SURE TO REDO THE CODING REGION VALUES, ADD IF STATEMENTS TO CHANGE WHAT GETS DISPLAYED BASED ON THE DROP DOWN MENU
    # volcano_fig, qq_fig, fig_mu, fig_sigma, dnds_fig, table_fig, text_special = generate_coding_region_report(dig_df, mutation_type, burden_type, display_labels_key, p_val_type)

    dig_data_columns = []
    
    for clm_nm in DIG_REPORT_COLUMN_NAMES:
        dig_data_clm_dict = {}
        dig_data_clm_dict["name"] = clm_nm
        dig_data_clm_dict["id"] = clm_nm

        if clm_nm == 'FDR' or "PVAL" in clm_nm:
            # gets the column name corresponding to the mutation type and the scatterpoint type
            if combined_burden_type[burden_type] == "":
                dig_data_clm_dict["id"] = clm_nm + "_"+ combined_mutation_type[mutation_type] + "_" + scatterpoint_type[p_val_type]
            else:
                dig_data_clm_dict["id"] = clm_nm + "_"+ combined_mutation_type[mutation_type] + "_" + combined_burden_type[burden_type] + "_" + scatterpoint_type[p_val_type]

        dig_data_columns.append(dig_data_clm_dict)
    
    for column_dict in dig_data_columns:
        # gets the dig data column name
        column = column_dict["id"]

        # skip the rank column
        if column == 'RANK':
            continue

        format='{:.3E}'

        # rounds all the values in the FDR and PVAL columns to 4 significant digits
        if 'FDR' in column or 'PVAL' in column:
            dig_df[column] = [format.format(value) for value in dig_df[column]]
        
    # get the coding region working plots working!!
    # get the display bounds selection tool working 

    # ONLY GETTING THE FIRST 500 ROWS OF DATA TO DISPLAY IN THE TABLE
    # REMOVE THE DEBUGGING LATER!!!   

    # sort the table with respect to 'PVAL' column, smallest(1) -> largest(nth)
    dig_df = dig_df[:100]

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
    dig_report_layout = [
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
        html.Div(children=[])
        ] + gen_dig_app_combined_component_layout()
    
    return dig_report_layout

def gen_dig_app_combined_component_layout(radio_label='Combined'):
    """ 
    
    """
    
    if radio_label == 'Combined':
        combined_layout = [
                # displays the interactive component to filter the samples displays based on their purity values
                html.Div([
                    dbc.Label(children="Debugging Stuff!!!", id="debugging"),    
                ]),

                # Plotly Figure for the DIG Report
                html.Div([
    
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
                            id='dig-report-combined-table',
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
    else:
        combined_layout = []

    return combined_layout

def gen_dig_app_coding_region_component_layout(radio_label='Coding regions'):
    """
    
    """
    if radio_label == 'Coding regions':
   
        coding_region_layout = [
                # displays the interactive component to filter the samples displays based on their purity values
                html.Div([
                    dbc.Label(children="Debugging Stuff!!!", id="debugging-coding"),    
                ]),

                # Plotly Figure for the DIG Report
                html.Div([

                    

                    dbc.Row([# insert the dropdown menus as columns inside this list for dbc.Row
                        dbc.Col([
                            # dropdown for selecting a mutation type
                            dbc.Label("Select Mutation Type"),
                            dcc.Dropdown(
                            id='dig-coding-mutation-dropdown',
                            options=[
                                {'label': 'Indels + Nonsynonymous SNVs', 'value': 'indels_nonsynonymous_snvs'},
                                {'label': 'Indels', 'value':'indels'},
                                {'label': 'Nonsynonymous SNVs', 'value':'nonsynonymous_snvs'},
                                {'label': 'Missense SNVs', 'value': 'missense_snvs'},
                                {'label': 'Nonsense SNVs', 'value': 'nonsense_snvs'},
                                {'label': 'Truncating SNVs', 'value': 'truncating_snvs'},
                                {'label': 'Splice Site SNVs', 'value': 'splice_site_snvs'},
                                {'label': 'Synonymous SNVs', 'value': 'synonymous_snvs'}
                            ],
                            value='',
                            ),
                        ]),
                        dbc.Col([
                            # dropdown for selecting burden type
                            dbc.Label("Select Burden Type"),
                            dcc.Dropdown(
                            id='dig-coding-burden-dropdown',
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
                            id='dig-coding-p-value-dropdown',
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
                                # makes a toggle component
                                daq.BooleanSwitch(
                                id='display-coding-bounds-toggle-switch',
                                label='Display Bounds',
                                on=False),
                            ]),
                            dbc.Col([
                                # makes a toggle component
                                daq.BooleanSwitch(
                                id='display-coding-labels-toggle-switch',
                                label='Display Labels',
                                on=False),
                            ])
                        ])
                    ]),

                    html.Div([
                        dbc.Label(id="coding-special-text-output", children=""),
                    ]),

                    # PLOTS ABOVE THE CODING REGION TABLE
                    dbc.Row([
                        dbc.Col([
                            # creates the dig volcano plot
                            dcc.Graph(id='dig-coding-volcano-graph', figure={}),

                        ]),
                        dbc.Col([
                            # creates the dig QQ plot
                            dcc.Graph(id='dig-coding-qq-graph', figure={}),
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
                                    html.Label(children="Combined", id="dig-coding-report-type-label"), # initialize label to empty string
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
                                    "id": i} for i in CODING_REGION_DIG_REPORT_COLUMN_NAMES
                            ],
                            data=pd.DataFrame(columns=CODING_REGION_DIG_REPORT_COLUMN_NAMES).to_dict(
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

                    # PLOTS BELOW THE CODING REGION TABLE
                    dbc.Row([
                        dbc.Col([
                            # creates the dig fig mu plot
                            dcc.Graph(id='dig-coding-fig-mu-graph', figure={}),
                        ]),
                        dbc.Col([
                            # creates the dig fig sigma plot
                            dcc.Graph(id='dig-coding-fig-sigma-graph', figure={}),
                        ]),
                        dbc.Col([
                            # creates the dig dnds plot
                            dcc.Graph(id='dig-coding-dnds-graph', figure={}),
                        ])
                    ])
                ], 
            )
        ]
    else:
        coding_region_layout = []

    return coding_region_layout

def gen_combined_dig_report_app_component():
    """
    
    """
    
    return AppComponent(
        name='DIG Component',
        layout=gen_dig_app_combined_component_layout(),
        new_data_callback=gen_dig_app_component_data_internal_callback,
        internal_callback=gen_dig_app_component_data_external_callback,
        callback_input=[
            Input('dig-report-type-label', 'children'),
            Input('dig-report-type-radioitems', 'value'),
            Input('dig-mutation-dropdown', 'value'),
            Input('dig-burden-dropdown', 'value'),
            Input('dig-p-value-dropdown', 'value'),
            Input('display-bounds-toggle-switch', 'on'),
            Input('display-labels-toggle-switch', 'on')
        ],

        callback_output=[
            Output('dig-report-combined-table', 'data'),
            Output('dig-report-type-label', 'children'),
            Output('dig-qq-graph', 'figure'),
            Output('dig-report-combined-table', 'columns'),
            Output('dig-mutation-dropdown', 'value'),
            Output('dig-burden-dropdown', 'value'),
            Output('dig-p-value-dropdown', 'value'),
            Output('special-text-output', 'children'),
            Output('debugging', 'children'),
        ],
    )
