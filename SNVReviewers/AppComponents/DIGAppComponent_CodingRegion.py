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

from SNVReviewers.AppComponents.utils import generate_coding_region_report, coding_region_mutation_type, coding_region_burden_type, scatterpoint_type

DIG_CODING_REGION_REPORT_COLUMN_NAMES = ["RANK", "GENE", 'CHROM', 'LENGTH', "FDR", "PVAL", "OBS", 
                                         'EXP', 'MU', 'SIGMA', 'dNdS_OBS', 'dNdS_EXP', 'FLAG', 'CGC', 'PANCAN']
                    
CODING_REGION_MUT_DROPDOWN = [
                        {'label': 'Indels + Nonsynonymous SNVs', 'value': 'indels_nonsynonymous_snvs'},
                        {'label': 'Indels', 'value': 'indels'},
                        {'label': 'Nonsynonymous + SNVs', 'value': 'nonsynonymous_snvs'},
                        {'label': 'Missense SNVs', 'value': 'missense_snvs'},
                        {'label': 'Nonsense SNVs', 'value': 'nonsense_snvs'},
                        {'label': 'Truncating SNVs', 'value': 'truncating_snvs'},
                        {'label': 'Splice Site SNVs', 'value': 'splice_site_snvs'},
                        {'label': 'Synonymous + SNVs', 'value': 'synonymous_snvs'}
                    ]
CODING_REGION_BUR_DROPDOWN = [
                        {'label': 'Total', 'value': 'total'},
                        {'label': 'Sample-wise', 'value': 'sample_wise'}
                    ]
CODING_REGION_SCATTER_DROPDOWN = [
                        {'label': 'Uniform P-mid', 'value': 'uniform_p_mid'},
                        {'label': 'P-mid', 'value': 'p_mid'},
                    ]


DIG_REPORT_VALUES = ["Combined", "Coding regions", "Promoter regions", "5-prime UTRs", "3-prime UTRs"]
DIG_DATAFRAME_IDX = 0
DND_DATAFRAME_IDX = 1
MUTSIG_DATAFRAME_IDX = 2
SNV_DATA_COLUMN_NAME = "snv_data"
DIG_LABEL_IDX = 1

def gen_dig_coding_region_app_component_data_internal_callback(
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

    # check if mutation_type is an old dropdown menu option

    if mutation_type == "":
        mutation_type = 'indels_nonsynonymous_snvs'
        burden_type = 'total'
        p_val_type = 'uniform_p_mid'

    dig_df = data.df[SNV_DATA_COLUMN_NAME][0][DIG_DATAFRAME_IDX].copy() # gets the dig report data for a specific cohort

    # sort the table with respect to 'PVAL' column, smallest(1) -> largest(nth)
    if dig_label == "Coding regions":
        dig_df = dig_df.sort_values(by='PVAL'+ "_"+ coding_region_mutation_type[mutation_type] + "_" + coding_region_burden_type[burden_type]) #+ "_" + scatterpoint_type[p_val_type])
    
    dig_df['RANK'] = np.array([i+1 for i in range(len(dig_df))])

    debugging_component = ""
    display_bounds = display_toggle_value # whether to display the bounds
    
    dir_output = "./example_notebooks/data/dig_data"
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

    # MAKE SURE TO REDO THE CODING REGION VALUES, ADD IF STATEMENTS TO CHANGE WHAT GETS DISPLAYED BASED ON THE DROP DOWN MENU
    df_kept, volcano_fig, qq_fig, fig_mu, fig_sigma, dnds_fig, table_fig, text_special = generate_coding_region_report(dig_df, mutation_type, burden_type, display_labels_key, p_val_type)
    df_kept_new = df_kept[['dNdS_OBS', 'dNdS_EXP', 'GENE']]
    dig_data_columns = []
    
    for clm_nm in DIG_CODING_REGION_REPORT_COLUMN_NAMES:
        dig_data_clm_dict = {}
        dig_data_clm_dict["name"] = clm_nm
        dig_data_clm_dict["id"] = clm_nm

        if "PVAL" in clm_nm:
            # gets the column name corresponding to the mutation type and the scatterpoint type
            dig_data_clm_dict["id"] = clm_nm + "_"+ coding_region_mutation_type[mutation_type] + "_" + coding_region_burden_type[burden_type] #+ "_" +scatterpoint_type[p_val_type]

        elif 'FDR' in clm_nm:
            dig_data_clm_dict["id"] = clm_nm + "_"+ coding_region_mutation_type[mutation_type] + "_" + scatterpoint_type[p_val_type]
        
        elif 'LENGTH' in clm_nm:
            dig_data_clm_dict['id'] = 'GENE_' +  clm_nm

        elif 'OBS' in clm_nm and coding_region_mutation_type[mutation_type] != 'MUT' and coding_region_mutation_type[mutation_type] != 'dNdS':
            dig_data_clm_dict["id"] = clm_nm + "_"+ coding_region_mutation_type[mutation_type]

        dig_data_columns.append(dig_data_clm_dict)
    
    for column_dict in dig_data_columns:
        # gets the dig data column name
        column = column_dict["id"]

        # skip the rank column
        if column == 'RANK':
            continue

        format='{:.3E}'

        # rounds all the values in the FDR and PVAL columns to 4 significant digits
        if 'FDR' in column or 'PVAL' in column or 'MU' in column or 'SIGMA' in column: #or 'dNdS_OBS' in column or 'dNdS_EXP' in column:
            dig_df[column] = [format.format(value) for value in dig_df[column]]
        
    # get the coding region working plots working!!
    # get the display bounds selection tool working 

    # ONLY GETTING THE FIRST 100 ROWS OF DATA TO DISPLAY IN THE TABLE
    # REMOVE THE DEBUGGING LATER!!!   
    dig_df = dig_df[:100]

    return [
            dig_df.to_dict('records'),
            dig_type_selection,
            volcano_fig,
            qq_fig,
            fig_mu, 
            fig_sigma,
            dnds_fig,
            dig_data_columns,
            mutation_type,
            burden_type,
            p_val_type,
            CODING_REGION_MUT_DROPDOWN,
            CODING_REGION_BUR_DROPDOWN,
            CODING_REGION_SCATTER_DROPDOWN, 

            # changing the plot size and the visibility of the plots
            {'width':'550px'},
            # {'display':'inline-block'}, # displays the volcano plot
            {'display':'inline-block', 'width':'550px'}, # displays the volcano plot
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
    output = gen_dig_coding_region_app_component_data_internal_callback(
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

def gen_dig_app_coding_region_component_layout():
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
                    value="Coding regions",
                    id="dig-report-type-radioitems",
                ),

                dbc.Row([# insert the dropdown menus as columns inside this list for dbc.Row
                    dbc.Col([
                        # dropdown for selecting a mutation type
                        dbc.Label("Select Mutation Type"),
                        dcc.Dropdown(
                        id='dig-mutation-dropdown',
                        options=CODING_REGION_MUT_DROPDOWN,
                        value='',
                        ),
                    ]),
                    dbc.Col([
                        # dropdown for selecting burden type
                        dbc.Label("Select Burden Type"),
                        dcc.Dropdown(
                        id='dig-burden-dropdown',
                        options=CODING_REGION_BUR_DROPDOWN,
                        value=''
                        ),
                    ]),
                    dbc.Col([
                        # dropdown for selecting burden type
                        dbc.Label("P-value Type"),
                        dcc.Dropdown(
                        id='dig-p-value-dropdown',
                        options=CODING_REGION_SCATTER_DROPDOWN,
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
                # Warning for doing specific dropdown menu combinations
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
                        # creates the dig volcano plot
                        dcc.Graph(id='dig-volcano-graph', 
                                  figure={},
                                  ),
                    ])
                ])
            ]),

            html.Div(
                [
                    # displays the type of dig report you want displayed
                    dbc.Row([
                        html.Div(
                            [
                                dbc.Label("Dig Coding Region Report Table: "),
                                html.Label(children="Coding Region", id="dig-report-type-label"), # initialize label to empty string
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
                                "id": i} for i in DIG_CODING_REGION_REPORT_COLUMN_NAMES
                        ],
                        data=pd.DataFrame(columns=DIG_CODING_REGION_REPORT_COLUMN_NAMES).to_dict(
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
                        dcc.Graph(id='dig-fig-mu-graph', figure={}),
                    ]),
                    dbc.Col([
                        # creates the dig fig sigma plot
                        dcc.Graph(id='dig-fig-sigma-graph', figure={}),
                    ]),
                    dbc.Col([
                        # creates the dig dnds fig plot
                        dcc.Graph(id='dig-dnds-fig-graph', figure={}),
                    ])
                ])
            ], 
        )
    ]

def gen_dig_coding_region_app_component():
    """
    
    """
    
    return AppComponent(
        name='DIG Coding Region Component',
        layout=gen_dig_app_coding_region_component_layout(),
        new_data_callback=gen_dig_coding_region_app_component_data_internal_callback,
        internal_callback=gen_dig_app_component_data_external_callback,
        callback_input=[
            # Input('dig-report-coding-region-type-label', 'children'),
            # Input('dig-report-coding-region-type-radioitems', 'value'),
            # Input('dig-coding-region-mutation-dropdown', 'value'),
            # Input('dig-coding-region-burden-dropdown', 'value'),
            # Input('dig-coding-region-p-value-dropdown', 'value'),
            # Input('display-bounds-coding-region-toggle-switch', 'on'),
            # Input('display-labels-coding-region-toggle-switch', 'on')
            Input('dig-report-type-label', 'children'),
            Input('dig-report-type-radioitems', 'value'),
            Input('dig-mutation-dropdown', 'value'),
            Input('dig-burden-dropdown', 'value'),
            Input('dig-p-value-dropdown', 'value'),
            Input('display-bounds-toggle-switch', 'on'),
            Input('display-labels-toggle-switch', 'on')
        ],

        callback_output=[
            # Output('dig-report-coding-region-table', 'data'),
            # Output('dig-report-coding-region-type-label', 'children'),
            # Output('dig-volcano-coding-region-graph', 'figure'),
            # Output('dig-qq-coding-region-graph', 'figure'),
            # Output('dig-fig-mu-coding-region-graph', 'figure'),
            # Output('dig-fig-sigma-coding-region-graph', 'figure'),
            # Output('dig-dnds-fig-coding-region-graph', 'figure'),
            # Output('dig-report-coding-region-table', 'columns'),
            # Output('dig-coding-region-mutation-dropdown', 'value'),
            # Output('dig-coding-region-burden-dropdown', 'value'),
            # Output('dig-coding-region-p-value-dropdown', 'value'),
            # Output('coding-region-special-text-output', 'children'),
            # Output('coding-region-debugging', 'children'),
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
