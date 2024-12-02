from dash import dcc, html
from dash.dependencies import Input, Output
import dash_bootstrap_components as dbc
import plotly.graph_objects as go

from SNVReviewers.AppComponents.utils import gen_cnp_figure
from AnnoMate.ReviewDataApp import AppComponent
from AnnoMate.DataTypes.GenericData import GenericData

def gen_mutsig_app_component_data_callback(
    data: GenericData,
    idx,
):
    
    return [idx]

def gen_mutsig_app_component_layout():
    
    # table
    #
    return html.Div(id='custom-component')


def gen_custom_app_component():
    
    return AppComponent(
        name='DIG Component',
        layout=gen_mutsig_app_component_layout(),
        new_data_callback=gen_mutsig_app_component_data_callback,
        internal_callback=gen_mutsig_app_component_data_callback,
        callback_output=[
            Output('custom-component', 'children'),
        ],
    )
