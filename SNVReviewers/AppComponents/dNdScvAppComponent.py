from dash import dcc, html
import dash_bootstrap_components as dbc
from AnnoMate.ReviewDataApp import AppComponent
from AnnoMate.DataTypes.GenericData import GenericData

from dash import dcc, html, dash_table
from dash.dependencies import Input, Output
import dash_bootstrap_components as dbc

from AnnoMate.ReviewDataApp import AppComponent
from AnnoMate.DataTypes.GenericData import GenericData
from SNVReviewers.AppComponents.dNdScvComponentHelpers import gen_dndscv_results_app_component, gen_dndscv_comparison_app_component, gen_dndscv_summary_app_component
from SNVReviewers.AppComponents.DIGAppComponent import DND_PLOT_DATAFRAME_IDX, DND_MERGED_DATAFRAME_IDX, DND_GLOBAL_DATAFRAME_IDX, DND_COMPARISON_DATAFRAME_IDX
from SNVReviewers.AppComponents.DIGAppComponent import SNV_DATA_COLUMN_NAME, PAGE_SIZE
from SNVReviewers.AppComponents.dNdScvComponentHelpers import RESULTS_DROPDOWN_OPTIONS, COMPARISON_DROPDOWN_OPTIONS, SUMMARY_DROPDOWN_OPTIONS, DNDSCV_REPORT_COLUMN_NAMES, DNDSCV_COMPARISON_COLUMN_NAMES, DNDSCV_SUMMARY_COLUMN_NAMES

import pandas as pd
import numpy as np

DNDSCV_REPORT_VALUES = ["Results", "Comparison", "Summary"]

def gen_dNdScv_app_component_data_callback(
    data: GenericData,
    idx,
    dnd_radio_item_selection,
    dnds_dropdown_value,
):
    

    # WORK ON THE DIG REPORT lollipop plot 
    # https://bioconductor.org/packages/release/bioc/vignettes/trackViewer/inst/doc/lollipopPlot.html
    # make the lollipop plot only for the DIG Report
    # also use Liz's plot as a base for getting the lollipop plot
        # https://github.com/getzlab/Lollipop_plots/blob/main/play_with_lollipop_plotting.ipynb
            # not in Python

    # Get the lollipop plot working for DIG report -> coding region and intro region!!

    all_page_content = []
    results_dropdown_values = [dict_option['value'] for dict_option in RESULTS_DROPDOWN_OPTIONS]
    comparison_dropdown_values = [dict_option['value'] for dict_option in COMPARISON_DROPDOWN_OPTIONS]
    summary_dropdown_values = [dict_option['value'] for dict_option in SUMMARY_DROPDOWN_OPTIONS]

    dnd_df_plot = data.df[SNV_DATA_COLUMN_NAME][0][DND_PLOT_DATAFRAME_IDX]
    dnd_df_merged = data.df[SNV_DATA_COLUMN_NAME][0][DND_MERGED_DATAFRAME_IDX]
    dnd_df_global = data.df[SNV_DATA_COLUMN_NAME][0][DND_GLOBAL_DATAFRAME_IDX]
    dnd_df_comparison = data.df[SNV_DATA_COLUMN_NAME][0][DND_COMPARISON_DATAFRAME_IDX]
    
    if dnds_dropdown_value is None:
        dnds_dropdown_value = 'all'

    if dnd_radio_item_selection == "Results":
        if dnds_dropdown_value not in results_dropdown_values:
            dnds_dropdown_value = 'all'

        all_page_content = gen_dndscv_results_app_component(
                dnd_df_plot,
                dnd_df_merged,
                dnd_df_global,
                dnds_dropdown_value,
                dnd_radio_item_selection
        )
    
    elif dnd_radio_item_selection == "Comparison":

        if dnds_dropdown_value not in comparison_dropdown_values:
            dnds_dropdown_value = "MutSig2CV vs dNdScv"

        all_page_content = gen_dndscv_comparison_app_component(
                dnd_df_comparison,
                dnds_dropdown_value,
                dnd_radio_item_selection
        )
        
    elif dnd_radio_item_selection == "Summary":

        if dnds_dropdown_value not in summary_dropdown_values:
            dnds_dropdown_value = "MutSig2CV, dNdScv and DIG"

        all_page_content = gen_dndscv_summary_app_component(
                dnd_df_comparison,
                dnds_dropdown_value,
                dnd_radio_item_selection
        )

    return all_page_content

def gen_dNdScv_app_component_layout():
    
    return [   
            # REMOVE LATER!!!
            html.Div([
                dbc.Label(id="dnds-debugging", children=""),
            ]),

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
                        html.Div([
                            dbc.Label(id="dnds-dropdowm-label", children=""),
                        ]),
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
                    dbc.Col([
                        # creates the dig QQ plot
                        dcc.Graph(id='dnds-comparison-graph', 
                                  figure={},
                                  style={"display":"none"} # increases the size of the plot
                                ), 
                    ]),
                ])
            ]),

            html.Div([
                    dbc.Label("")
                ]),

            html.Div(
                [        
                # displays first table for the dNdScv comparison report
                html.Div(
                    children=[
                        dash_table.DataTable(
                            id='dnds-comparison-table1',
                            columns=[
                                {"name": i,
                                    "id": i} for i in DNDSCV_COMPARISON_COLUMN_NAMES
                            ],
                            data=pd.DataFrame(columns=DNDSCV_COMPARISON_COLUMN_NAMES).to_dict(
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
                            page_size=PAGE_SIZE,
                        
                            # changing the width of the data table to 
                            style_table={
                                'display':"none"
                                },
                            style_data_conditional=[]
                        ),
                    ]
                ),
                html.Div([
                    dbc.Label("")
                ]),

                html.Div(
                    children=[
                        dash_table.DataTable(
                        id='dnds-comparison-table2',
                        columns=[
                            {"name": i,
                                "id": i} for i in DNDSCV_COMPARISON_COLUMN_NAMES
                        ],
                        data=pd.DataFrame(columns=DNDSCV_COMPARISON_COLUMN_NAMES).to_dict(
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
                        page_size=PAGE_SIZE,
                    
                        # changing the width of the data table to 
                        style_table={
                            'display':"none"
                            },
                        style_data_conditional=[]
                        ),
                    ]
                ),

                html.Div([
                    dbc.Label("")
                ]),

                html.Div(
                    children=[
                        dash_table.DataTable(
                        id='dnds-comparison-table3',
                        columns=[
                            {"name": i,
                                "id": i} for i in DNDSCV_COMPARISON_COLUMN_NAMES
                        ],
                        data=pd.DataFrame(columns=DNDSCV_COMPARISON_COLUMN_NAMES).to_dict(
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
                        page_size=PAGE_SIZE,
                    
                        # changing the width of the data table to 
                        style_table={
                            'display':"none"
                            },
                        style_data_conditional=[]
                        ),
                    ]
                ),
                html.Div([
                    dbc.Label("")
                ]),
                # displays a table for the dig report
                html.Div(
                    children=[
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
                        page_size=PAGE_SIZE,
                    
                        # changing the width of the data table to 
                        style_table={
                            'width': '100%',  # Make the table width responsive
                            'maxWidth': '100%',  # Ensure it doesn’t go beyond the screen width
                            'overflowX': 'auto',  # Allow horizontal scroll if necessary
                            },
                        style_data_conditional=[]
                        ),
                    ]
                ),
                html.Div(
                    children=[
                        dash_table.DataTable(
                        id='dnds-summary-table',
                        columns=[
                            {"name": i,
                                "id": i} for i in DNDSCV_SUMMARY_COLUMN_NAMES
                        ],
                        data=pd.DataFrame(columns=DNDSCV_SUMMARY_COLUMN_NAMES).to_dict(
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
                        page_size=PAGE_SIZE,
                    
                        # changing the width of the data table to 
                        style_table={
                            'display':"none"
                            },
                        style_data_conditional=[]
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
            Input('dnds-report-type-radioitems', 'value'),
            Input('dnds-gene-dropdown', 'value')
        ],
        
        callback_output=[
            # dNdScv tables
            Output('dnds-report-table', 'data'),
            Output('dnds-comparison-table1', 'data'),
            Output('dnds-comparison-table2', 'data'),
            Output('dnds-comparison-table3', 'data'),
            # Output('dnds-summary-table', 'figure'),
            Output('dnds-summary-table', 'data'),

            # updates the columns for the comparison tables
            # Output('dnds-report-table', 'columns'),
            Output('dnds-comparison-table1', 'columns'),
            Output('dnds-comparison-table2', 'columns'),
            Output('dnds-comparison-table3', 'columns'),
            Output('dnds-summary-table', 'columns'),

            # dNdScv figures
            Output('dnds-qq-graph', 'figure'),
            Output('dnds-mutation-ratio-graph', 'figure'),
            Output('dnds-missense-graph', 'figure'),
            Output('dnds-truncating-graph', 'figure'),
            Output('dnds-comparison-graph', 'figure'),

            # dNdScv figure display
            Output('dnds-qq-graph', 'style'),
            Output('dnds-mutation-ratio-graph', 'style'),
            Output('dnds-missense-graph', 'style'),
            Output('dnds-truncating-graph', 'style'),
            Output('dnds-comparison-graph', 'style'),

            # dNdScv table styling for sizing, hiding, or displaying tables
            Output('dnds-report-table', 'style_table'),
            Output('dnds-comparison-table1', 'style_table'),
            Output('dnds-comparison-table2', 'style_table'),
            Output('dnds-comparison-table3', 'style_table'),
            Output('dnds-summary-table', 'style_table'),

            # dNdScv table conditional styling of data for significant genes
            Output('dnds-report-table', 'style_data_conditional'),
            Output('dnds-comparison-table1', 'style_data_conditional'),
            Output('dnds-comparison-table2', 'style_data_conditional'),
            Output('dnds-comparison-table3', 'style_data_conditional'),
            Output('dnds-summary-table', 'style_data_conditional'),

            Output('dnds-special-text-output', 'children'),
            Output('dnds-gene-dropdown', 'options'),
            Output('dnds-dropdowm-label', 'children'),
            Output('dnds-gene-dropdown', 'value'),

            Output('dnds-debugging', 'children'),
        ],
    )