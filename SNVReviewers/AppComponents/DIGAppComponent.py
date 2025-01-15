from dash import dcc, html
from dash.dependencies import Input, Output
import dash_bootstrap_components as dbc
import plotly.graph_objects as go

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
    p_val_type
):
    """
    
    """

    if mutation_type == "":
        mutation_type = 'indels_snvs'
        burden_type = 'total'
        p_val_type = 'uniform_p_mid'

    # add a radio item that allows the user to be able to select which type of dig table they want to view
        # coding
        # 3-utr
        # 5-utr
        # non coding
        # combined
    dig_df = data.df[SNV_DATA_COLUMN_NAME][0][DIG_DATAFRAME_IDX] # gets the dig report data for a specific cohort
    dnd_df = data.df[SNV_DATA_COLUMN_NAME][0][DND_DATAFRAME_IDX]
    mutsig_df = data.df[SNV_DATA_COLUMN_NAME][0][MUTSIG_DATAFRAME_IDX]
    debugging_component = "" #+ str(type(data.df["snv_data"]))
    dig_output_type = "Combined"
    display_bounds = True # whether to display the bounds
    
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

        debugging_component += "name: " + dig_data_clm_dict["name"] + " id: " + dig_data_clm_dict["id"] + "; "
        
        dig_data_columns.append(dig_data_clm_dict)
    

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
    p_val_type
):
    """
    """
    # add a radio item that allows the user to be able to select which type of dig table they want to view
        # coding
        # 3-utr
        # 5-utr
        # non coding
        # combined
    # print(type)
    output = gen_dig_app_component_data_internal_callback(
                data,
                idx,
                dig_label,
                dig_type_selection, 
                mutation_type,
                burden_type,
                p_val_type
            )
    
    # output[DIG_LABEL_IDX] = dig_label

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

                # dropdown for selecting a mutation type
                dbc.Label("Select Mutation Type"),
                dcc.Dropdown(
                id='dig-mutation-dropdown',
                options=[
                    {'label': 'Indels + SNVs', 'value': 'indels_snvs'},
                    {'label': 'Indels', 'value': 'indels'},
                    {'label': 'SNVs', 'value': 'snvs'}
                ],
                value=''
                ),


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

                # dropdown for selecting burden type
                dbc.Label("P-value Type"),
                dcc.Dropdown(
                id='dig-p-value-dropdown',
                options=[
                    {'label': 'Uniform P-mid', 'value': 'uniform_p_mid'},
                    {'label': 'P-mid', 'value': 'p_mid'},
                ],
                value=''),

                dbc.Label(id="special-text-output", children=""),

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
                        page_size=5),
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


# PVAL_coding_SNV_recalc	
    # PVAL_coding_SNV_unif	
    # PVAL_coding_SNV_lower	
    # PVAL_coding_SNV_upper	
    # PVAL_coding_SNV_SAMPLE_recalc	
    # PVAL_coding_SNV_SAMPLE_unif	
    # PVAL_coding_SNV_SAMPLE_lower	
    # PVAL_coding_SNV_SAMPLE_upper	
    # PVAL_coding_INDEL_recalc	
    # PVAL_coding_INDEL_unif	
    # PVAL_coding_INDEL_lower	
    # PVAL_coding_INDEL_upper	
    # PVAL_coding_MUT_recalc	
    # PVAL_coding_MUT_unif	
    # PVAL_coding_MUT_lower	
    # PVAL_coding_MUT_upper	
    # SIZE_promoter	
    # PVAL_promoter_SNV_recalc	
    # PVAL_promoter_SNV_unif	
    # PVAL_promoter_SNV_lower	
    # PVAL_promoter_SNV_upper	
    # PVAL_promoter_SNV_SAMPLE_recalc	
    # PVAL_promoter_SNV_SAMPLE_unif	
    # PVAL_promoter_SNV_SAMPLE_lower	
    # PVAL_promoter_SNV_SAMPLE_upper	
    # PVAL_promoter_INDEL_recalc	
    # PVAL_promoter_INDEL_unif	
    # PVAL_promoter_INDEL_lower	
    # PVAL_promoter_INDEL_upper	
    # PVAL_promoter_MUT_recalc	
    # PVAL_promoter_MUT_unif	
    # PVAL_promoter_MUT_lower	
    # PVAL_promoter_MUT_upper	
    # SIZE_5utr	
    # PVAL_5utr_SNV_recalc	
    # PVAL_5utr_SNV_unif	PVAL_5utr_SNV_lower	PVAL_5utr_SNV_upper	PVAL_5utr_SNV_SAMPLE_recalc	PVAL_5utr_SNV_SAMPLE_unif	
    # PVAL_5utr_SNV_SAMPLE_lower	PVAL_5utr_SNV_SAMPLE_upper	PVAL_5utr_INDEL_recalc	PVAL_5utr_INDEL_unif	
    # PVAL_5utr_INDEL_lower	PVAL_5utr_INDEL_upper	PVAL_5utr_MUT_recalc	PVAL_5utr_MUT_unif	PVAL_5utr_MUT_lower	
    # PVAL_5utr_MUT_upper	SIZE_3utr	PVAL_3utr_SNV_recalc	PVAL_3utr_SNV_unif	PVAL_3utr_SNV_lower	PVAL_3utr_SNV_upper	
    # PVAL_3utr_SNV_SAMPLE_recalc	PVAL_3utr_SNV_SAMPLE_unif	PVAL_3utr_SNV_SAMPLE_lower	PVAL_3utr_SNV_SAMPLE_upper	
    # PVAL_3utr_INDEL_recalc	PVAL_3utr_INDEL_unif	PVAL_3utr_INDEL_lower	PVAL_3utr_INDEL_upper	PVAL_3utr_MUT_recalc	
    # PVAL_3utr_MUT_unif	PVAL_3utr_MUT_lower	PVAL_3utr_MUT_upper	PVAL_SNV_recalc	FDR_SNV_recalc	PVAL_SNV_unif	FDR_SNV_unif	
    # PVAL_SNV_lower	FDR_SNV_lower	PVAL_SNV_upper	FDR_SNV_upper	PVAL_INDEL_recalc	FDR_INDEL_recalc	PVAL_INDEL_unif	
    # FDR_INDEL_unif	PVAL_INDEL_lower	FDR_INDEL_lower	PVAL_INDEL_upper	FDR_INDEL_upper	PVAL_MUT_recalc	FDR_MUT_recalc	
    # PVAL_MUT_unif	FDR_MUT_unif	PVAL_MUT_lower	FDR_MUT_lower	PVAL_MUT_upper	FDR_MUT_upper	PVAL_SNV_SAMPLE_recalc	
    # FDR_SNV_SAMPLE_recalc	PVAL_SNV_SAMPLE_unif	FDR_SNV_SAMPLE_unif	PVAL_SNV_SAMPLE_lower	FDR_SNV_SAMPLE_lower	
    # PVAL_SNV_SAMPLE_upper	FDR_SNV_SAMPLE_upper	CGC	PANCAN