from dash import dcc, html
from dash.dependencies import Input, Output
import dash_bootstrap_components as dbc
import plotly.graph_objects as go

from PurityReviewers.AppComponents.utils import gen_cnp_figure
from AnnoMate.ReviewDataApp import AppComponent
from AnnoMate.DataTypes.GenericData import GenericData

def gen_custom_app_component():
    
    return AppComponent(
        name='My Custom Component',
        layout=gen_custom_app_component_layout(),
        new_data_callback=custom_app_component_callback,
        callback_output=[
            Output('custom-component', 'children'),
        ],
    )

def gen_custom_app_component_layout():
    
    return html.Div(id='custom-component')

def custom_app_component_callback(
    data: GenericData,
    idx,
)
    return [idx]
