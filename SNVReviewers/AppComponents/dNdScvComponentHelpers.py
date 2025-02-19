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

# from SNVReviewers.AppComponents.utils import ngenes


RESULTS_DROPDOWN_VALUES = ["All", "First 30", "First 20", "First 10", "None"]
COMPARISON_DROPDOWN_VALUES = ["MutSig2 vs dNdScv", "MutSig2 vs DIG", "dNdScv vs DIG"]
SUMMARY_DROPDOWN_VALUES = ["MutSig2, dNdScv, and DIG", "MutSig2 and dNdScv only", "MutSig2 and DIG only",
                           "dNdScv and DIG only", "MutSig2 only", "dNdScv only", "DIG only"]

# dNdScv dataframe and plot generation
import argparse
import pandas as pd
import plotly.graph_objects as go
import numpy as np
import json

#
# Visualization Parameters

# minimum number of rows (genes) to display in the table
n_rows_min = 50
# horizontal buffer for the scatter plots
hor_buffer = 0.05
hor_buffer_global = 0.25
# vertical gap between gene labels and scatter points
y_gap_annot = 0.01
# buffer for the number of rows in the table
n_rows_buffer = 0.5
# maximum value along the axis of p-values
pval_max = 16

# properties of significant points
col_sig = 'rgba(255, 0, 0, 1)'
opac_sig = 0.7
# properties of near-significant points
col_nearsig = 'rgba(30, 144, 255, 1)'
opac_nearsig = 0.7
# properties of non-significant points
col_nonsig = 'rgba(0, 0, 0, 1)'
opac_nonsig = 0.5
# properties of thinner lines
thk_thin = 0.75
col_thin = 'gray'
typ_thin = 'dash'
# properties of thicker lines
thk_thick = 2.5
col_thick = 'gray'
typ_thick = 'dash'
# properties of error bars
wid_err = 1
wid_err_global = 2
thk_err = 0.5
thk_err_global = 1
opac_err = 0.3
# font size of axes labels
fs_axes = 14

# markers used to indicate the dominant test type in QQ plots
markers = {
    'mis': ['Missense', 'circle'],
    'trunc': ['Truncating', 'square'],
    'ind': ['InDel', 'star']
}
msize = 7.5

# dropdown options: number of genes to display in the scatter plots
ngenes = {
    'All': 'all',
    'First 30': 30,
    'First 20': 20,
    'First 10': 10,
    'None': None
}

# derived parameters
col_err_sig = ','.join(col_sig.split(',')[:-1]) + ', {})'.format(opac_err)
col_err_nearsig = ','.join(col_nearsig.split(',')[:-1]) + ', {})'.format(opac_err)
col_err_nonsig = ','.join(col_nonsig.split(',')[:-1]) + ', {})'.format(opac_err)


def reformat_numbers(x, format='{:.2E}'):
    """
    Reformat numbers in an array to a specific format
    """
    return [format.format(n) for n in x]


def get_capped_labels(x_capped, y_capped, labels_capped):
    labels_capped = [(
        f"({x_capped[i]:.2f}, {y_capped[i]:.2f})<br>{l}")
        for i, l in enumerate(labels_capped)]
    return np.array(labels_capped)


def prepare_scatter_data(labels, x, y, y_lower, y_upper, ind, col_err, thk_err, wid_err, capped=False):
    """
    Prepare scatter plot data and error bars for given indices.

    Parameters:
    labels (np.ndarray): Array of labels.
    x (np.ndarray): Array of x values.
    y (np.ndarray): Array of y values.
    y_lower (np.ndarray): Array of lower y values.
    y_upper (np.ndarray): Array of upper y values.
    ind (np.ndarray): Boolean array indicating the indices to use.
    col_err (str): Color for the error bars.
    thk_err (float): Thickness of the error bars.
    wid_err (float): Width of the error bars.
    labels_selected (np.ndarray, optional): Array of selected labels. Defaults to None.

    Returns:
    dict: Dictionary containing scatter plot data and error bars.
    """
    x_selected = x[ind]
    y_selected = y[ind]
    y_selected_lower = y_lower[ind]
    y_selected_upper = y_upper[ind]
    y_selected_upper_err = y_selected_upper - y_selected
    y_selected_lower_err = y_selected - y_selected_lower
    if capped:
        labels_selected = labels[ind]
        labels_selected = np.array([(
            f"({x_selected[i]:.2f}, {y_selected[i]:.2f} +{y_selected_upper_err[i]:.2f} / -{y_selected_lower_err[i]:.2f})<br>{l}")
            for i, l in enumerate(labels_selected)])
        x[ind] = pval_max
        x_selected = np.array([pval_max] * ind.sum())
    else:
        labels_selected = labels[ind]
    dict_erry = dict(
        type='data',
        symmetric=False,
        array=np.round(y_selected_upper_err, 2).tolist(),
        arrayminus=np.round(y_selected_lower_err, 2).tolist(),
        thickness=thk_err,
        width=wid_err,
        color=col_err
    )
    return {
        'labels': labels_selected,
        'x': x_selected,
        'y': y_selected,
        'dict_erry': dict_erry
    }


def plot_dnds(df, mut_typ, mut_name, alp, alp_nearsig):
    """
    Generate a scatter plot with FDR values on the horizontal and dNdS values on the vertical axis.

    Parameters:
    df (pd.DataFrame): DataFrame containing dNdS values and confidence intervals.
    mut_typ (str): Mutation type (e.g., 'mis' for missense, 'tru' for truncating).
    mut_name (str): Mutation name for labeling the plot.
    alp (float): Significance level for determining significant genes.
    alp_nearsig (float): Near-significance threshold.

    Returns:
    plotly.graph_objects.Figure: Plotly figure object containing the scatter plot.
    """
    labels = df.index.copy().to_numpy()
    x = -np.log10(df.qglobal_cv.copy().to_numpy())
    y = df[mut_typ + '_mle'].copy().to_numpy()
    y_lower = df[mut_typ + '_low'].copy().to_numpy()
    y_upper = df[mut_typ + '_high'].copy().to_numpy()
    # indicator of significant genes
    ind_sig = x >= -np.log10(alp)
    # indicator of capped FDR values
    ind_capped = x > pval_max
    # indicator of significant genes that are not capped
    ind_ncapped = np.logical_and(ind_sig, ~ind_capped)
    # indicator of significant or near-significant genes
    ind_hits = x >= -np.log10(alp_nearsig)
    # indicator of near-significant genes
    ind_nearsig = np.logical_and(~ind_sig, ind_hits)
    # indicator of non-significant genes
    ind_nsig = ~ind_sig

    # scatter points for non-significant genes
    nsig_data = prepare_scatter_data(labels, x, y, y_lower, y_upper, ind_nsig, col_err_nonsig, thk_err, wid_err)
    # scatter points for nearly-significant genes
    nearsig_data = prepare_scatter_data(labels, x, y, y_lower, y_upper, ind_nearsig, col_err_nearsig, thk_err, wid_err)
    # scatter points for significant genes with non-capped FDR values
    ncapped_data = prepare_scatter_data(labels, x, y, y_lower, y_upper, ind_ncapped, col_err_sig, thk_err, wid_err)
    # scatter points for significant genes with capped FDR values
    capped_data = prepare_scatter_data(labels, x, y, y_lower, y_upper, ind_capped, col_err_sig, thk_err, wid_err, capped=True)

    plotly_fig = go.Figure()

    scatter_data = [
        (nsig_data['x'].tolist(), nsig_data['y'].tolist(), nsig_data['dict_erry'], col_nonsig, opac_nonsig,
         nsig_data['labels'].tolist(), 'Non-significant'),
        (nearsig_data['x'].tolist(), nearsig_data['y'].tolist(), nearsig_data['dict_erry'], col_nearsig, opac_nearsig,
         nearsig_data['labels'].tolist(), 'Nearly-significant'),
        (ncapped_data['x'].tolist(), ncapped_data['y'].tolist(), ncapped_data['dict_erry'], col_sig, opac_sig,
         ncapped_data['labels'].tolist(), 'Significant'),
        (capped_data['x'].tolist(), capped_data['y'].tolist(), capped_data['dict_erry'], col_sig, opac_sig,
         capped_data['labels'], 'Significant')
    ]

    for i, (x_data, y_data, error_y, color, opacity, text, name) in enumerate(scatter_data):
        if i == len(scatter_data)-1:
            hoverinfo = 'name+text'
        else:
            hoverinfo = 'all'
        plotly_fig.add_trace(
            go.Scatter(
                x=x_data,
                y=y_data,
                error_y=error_y,
                mode='markers',
                marker=dict(color=color, opacity=opacity),
                text=text,
                name=name,
                showlegend=False,
                hoverinfo=hoverinfo,
                xhoverformat='.2f',
                yhoverformat='.2f'
            )
        )

    # define boundaries
    xlim = [0, max(x) * (1 + hor_buffer)]
    ylim = [0, max(max(y_upper[ind_sig]), max(y)) * (1 + hor_buffer)]
    # plot cap line
    plotly_fig.add_trace(
        go.Scatter(
            x=[pval_max] * 2,
            y=ylim,
            mode='lines',
            line=dict(dash=typ_thin, color=col_thin, width=thk_thin),
            showlegend=False,
            hoverinfo='skip'
        )
    )
    # format and add lables
    plotly_fig.update_layout(
        title=dict(
            text='dNdS' + ' Ratio of ' + (mut_name + ' mutations:').title(),
            font=dict(size=20)  # title font size
        ),
        xaxis_title='-Log10(FDR)',
        yaxis_title='dNdS ratio',
        xaxis=dict(
            range=xlim
        ),
        yaxis=dict(
            range=ylim
        ),
        template='plotly_white'
    )
    return plotly_fig


def plot_qq(x, y, labels, col, opac, mark, name, fig, hi='x+y+text+name'):
    fig.add_trace(
        go.Scatter(
            x=x.tolist(),
            y=y.tolist(),
            mode='markers',
            marker=dict(color=col, opacity=opac, symbol=mark, size=msize),
            text=labels.tolist(),
            name=name,
            hoverinfo=hi,
            showlegend=False,
            xhoverformat='.2f',
            yhoverformat='.2f'
        )
    )
    return fig


def generate_report(
        path_dnds_out,
        path_dnds_ci,
        path_dnds_global,
        path_cgc_list,
        path_pancan_list,
        dir_output,
        prefix_output=None,
        alp=0.1,
        alp_nearsig=0.25
):
    """
    Generate a dNdScv report for coding regions.

    Parameters:
    path_dnds_out (str): Path to the dNdScv output file containing p-values.
    path_dnds_ci (str): Path to the dNdScv confidence intervals file.
    path_dnds_global (str): Path to the dNdScv global values file.
    path_cgc_list (str): Path to the list of CGC genes.
    path_pancan_list (str): Path to the list of PanCanAtlas genes.
    dir_output (str): Output directory.
    prefix_output (str, optional): Prefix for the output file. Defaults to None.
    alp (float, optional): Significance level. Defaults to 0.1.
    alp_nearsig (float, optional): Near-significance threshold (default: 0.25).
    Returns:
    None
    """
    # lists of known driver genes
    cgc_list = pd.read_csv(path_cgc_list, sep='\t').to_numpy().flatten()
    pancan_list = pd.read_csv(path_pancan_list, sep='\t').to_numpy().flatten()
    # output from dNdScv containing p-values
    df_out = pd.read_csv(path_dnds_out, sep='\t').set_index('gene_name')
    # output from dNdScv containing dNdS values and confidence intervals
    df_ci = pd.read_csv(path_dnds_ci, sep='\t').set_index('gene')
    # merging the two dataframes
    df_merged = pd.concat([df_out, df_ci], axis=1)
    # global dNdS values
    df_global = pd.read_csv(path_dnds_global, sep='\t')
    # dataframe that will be plotted in a table format
    df_plot = df_merged.copy().reset_index().rename(columns={'index': 'GENE'})
    df_plot['RANK'] = df_plot.index + 1
    # adding indicator of genes being part of the CGC or PanCan list
    df_plot['CGC'] = df_plot.GENE.isin(cgc_list)
    df_plot['PANCAN'] = df_plot.GENE.isin(pancan_list)
    # renaming columns
    dict_n = {c: c.upper() for c in ['n_syn', 'n_mis', 'n_non', 'n_spl', 'n_ind']}
    dict_w = {c: 'dNdS_' + c.split('_')[0][1:].upper() for c in ['wmis_cv', 'wnon_cv', 'wspl_cv', 'wind_cv']}
    dict_pvals = {c: 'PVAL_' + c.split('_')[0][1:].upper() for c in ['pmis_cv', 'ptrunc_cv', 'pind_cv']}
    df_plot = df_plot.rename(columns=dict_n)
    df_plot = df_plot.rename(columns=dict_w)
    df_plot = df_plot.rename(columns=dict_pvals)
    df_plot = df_plot.rename(columns={'pglobal_cv': 'PVAL', 'qglobal_cv': 'FDR'})
    cols_kept = (['RANK', 'GENE'] + list(dict_n.values()) + list(dict_w.values()) + list(dict_pvals.values()) +
                 ['PVAL', 'FDR', 'CGC', 'PANCAN'])
    # cols_kept = (['RANK', 'GENE'] + list(dict_n.values()) + list(dict_w.values()) + ['PVAL', 'FDR', 'CGC', 'PANCAN'])

    #
    # HTML Template

    html_content = """
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>dNdScv Results</title>
        <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
        <style>
            .container {{
                display: flex;
            }}
            .plot {{
                margin: 10px;
            }}
            #qq-plot {{
                width: 100%;
            }}
            #table-plot {{
                width: 100%;
            }}
            #dnds-global {{
                width: 33%;
            }}
            #dnds-mis {{
                width: 33%;
            }}
            #dnds-tru {{
                width: 33%;
            }}
            .figures-container {{
                display: flex;
                justify-content: space-between;
            }}
            .figure-column {{
                width: 33%;
            }}
        </style>
    </head>
    <body>
        <h1>dNdScv Results</h1>
    
        <label for="ngenes">Number of Significant Gene Labels to Display:</label>
        <select id="ngenes" onchange="updatePlot()">
            {ngenes_options}
        </select>
    
        <div class="container">
            <div id="qq-plot" class="plot"></div>
        </div>
        <div class="table-plot">{fig_table_html}</div>
        <div class="figures-container">
            <div class="figure-column">{fig_dnds_global_html}</div>
            <div class="figure-column">{fig_dnds_mis_html}</div>
            <div class="figure-column">{fig_dnds_tru_html}</div>
        </div>
    
        <script>
            var plotData = {plot_data};
    
            function updatePlot() {{
                var ngenesTypeKey = document.getElementById("ngenes").value;
    
                var data = plotData[ngenesTypeKey];
    
                // Update Q-Q Plot
                var qqData = data.qq;
                Plotly.react('qq-plot', qqData);
            }}
    
            // Initial plot
            updatePlot();
        </script>
    </body>
    </html>
    """

    # generate the dropdown options
    ngenes_options = "\n".join([f'<option value="{key}">{key}</option>' for key in ngenes.keys()])

    tests_pvals = ['p' + mut_typ + '_cv' for mut_typ in markers.keys()]
    # prepare plot data for all dropdown options
    plot_data = {}
    for ngenes_key, ngenes_val in ngenes.items():

        #
        # Q-Q Plot

        # scatter plots
        labels = df_merged.index.copy().to_numpy()
        pvals = df_merged.pglobal_cv.copy().to_numpy()
        x = -np.log10(np.arange(1, len(pvals) + 1) / (len(pvals) + 1))
        y = -np.log10(pvals)
        logfdr = -np.log10(df_merged.qglobal_cv.copy().to_numpy())
        # indicator of "dominant" test (with minimal FDR) for each gene
        test_dom = df_merged[tests_pvals].idxmin(axis=1).to_numpy()
        # indicator of genes whose p-value exceeds the maximum value
        ind_capped = y > pval_max
        # indicator of significant genes
        ind_sig = logfdr >= -np.log10(alp)
        # indicator of significant or near-significant genes
        ind_hits = logfdr >= -np.log10(alp_nearsig)
        # indicator of near-significant genes
        ind_nearsig = np.logical_and(~ind_sig, ind_hits)
        # indicator of non-significant genes
        ind_nsig = ~ind_sig
        # indicator of hits that are not capped
        ind_ncapped = np.logical_and(ind_sig, ~ind_capped)
        # indicator of hits that are not capped
        ind_ncapped_hits = np.logical_and(ind_hits, ~ind_capped)
        ylim_upper = min(np.max(y), pval_max) * (1 + hor_buffer)

        # Scatter plot
        qq_fig = go.Figure()
        for i in range(len(tests_pvals)):
            test_i = tests_pvals[i][1:-3]
            ind = test_dom == tests_pvals[i]
            i_nsig = np.logical_and(ind, ind_nsig)
            i_nearsig = np.logical_and(ind, ind_nearsig)
            i_ncapped = np.logical_and(ind, ind_ncapped)
            i_capped = np.logical_and(ind, ind_capped)
            qq_fig = plot_qq(x[i_nsig], y[i_nsig], labels[i_nsig], col_nonsig, opac_nonsig, markers[test_i][1],
                             'Non-significant', qq_fig)
            qq_fig = plot_qq(x[i_nearsig], y[i_nearsig], labels[i_nearsig], col_nearsig, opac_nearsig,
                             markers[test_i][1],
                             'Near-significance', qq_fig)
            qq_fig = plot_qq(x[i_ncapped], y[i_ncapped], labels[i_ncapped], col_sig, opac_sig, markers[test_i][1],
                             'Significant', qq_fig)
            xi_capped = x[i_capped]
            qq_fig = plot_qq(xi_capped, np.array([pval_max] * sum(i_capped)),
                             get_capped_labels(xi_capped, y[i_capped], labels[i_capped]), col_sig, opac_sig,
                             markers[test_i][1], 'Significant', qq_fig, hi='name+text')

        # Dummy points for legend
        tests_dom, test_counts = np.unique(test_dom, return_counts=True)
        for test_pval_i in tests_dom[np.argsort(test_counts)[::-1]]:
            test_i = test_pval_i[1:-3]
            qq_fig.add_trace(
                go.Scatter(
                    y=[None],
                    mode='markers',
                    marker=dict(
                        color='white',
                        symbol=markers[test_i][1],
                        size=msize * 1.25,
                        line=dict(color='black', width=2)
                    ),
                    name=markers[test_i][0]
                )
            )

        # Add rotated text labels as annotations
        x_capped = x[ind_capped]
        y_capped = y[ind_capped]
        labels_capped = get_capped_labels(x_capped, y_capped, labels[ind_capped])
        if ngenes_val is not None:
            if ngenes_val == 'all':
                x_annot = [x_capped, x[ind_ncapped_hits]]
                y_annot = [[pval_max] * sum(ind_capped), y[ind_ncapped_hits]]
                labels_annot = [labels_capped, labels[ind_ncapped_hits]]
            else:
                if ngenes_val <= len(x_capped):
                    x_annot = [x_capped[:ngenes_val], []]
                    y_annot = [[pval_max] * ngenes_val, []]
                    labels_annot = [labels_capped[:ngenes_val], []]
                else:
                    end_ncapped = ngenes_val - len(x_capped)
                    x_annot = [x_capped, x[ind_ncapped_hits][:end_ncapped]]
                    y_annot = [[pval_max] * sum(ind_capped), y[ind_ncapped_hits][:end_ncapped]]
                    labels_annot = [labels_capped, labels[ind_ncapped_hits][:end_ncapped]]
            count = 0
            for (xa, ya, laba) in zip(x_annot, y_annot, labels_annot):
                for (xi, yi, label) in zip(xa, ya, laba):
                    if ind_nearsig[count]:
                        coli = col_nearsig
                    else:
                        coli = col_sig
                    qq_fig.add_annotation(
                        x=xi,
                        y=yi - ylim_upper * y_gap_annot,
                        text=label.split('<br>')[-1],
                        showarrow=False,
                        font=dict(color=coli),
                        textangle=-90,
                        xanchor="center",
                        yanchor="top"
                    )
                    count += 1

        # line plots
        qq_fig.add_trace(
            go.Scatter(
                x=[0, np.max(x) * (1 + hor_buffer)],
                y=[0, np.max(x) * (1 + hor_buffer)],
                mode='lines',
                line=dict(dash=typ_thick, color=col_thick, width=thk_thick),
                showlegend=False,
                hoverinfo='skip'
            )
        )
        qq_fig.add_trace(
            go.Scatter(
                x=[0, np.max(x) * (1 + hor_buffer)],
                y=[pval_max] * 2,
                mode='lines',
                line=dict(dash=typ_thin, color=col_thin, width=thk_thin),
                showlegend=False,
                hoverinfo='skip'
            )
        )
        # figure formatting
        qq_fig.update_layout(
            title='QQ-Plot of P-values:',
            xaxis_title='Expected -Log10(P-value)',
            yaxis_title='Observed -Log10(P-value)',
            xaxis=dict(range=[0, np.max(x) * (1 + hor_buffer)]),
            yaxis=dict(range=[0, ylim_upper]),
            template='plotly_white',
            legend=dict(
                title=dict(
                    text="Dominant Mutation Type:",
                ),
                indentation=10
            )
        )
        # save figures as separate data
        plot_data[f"{ngenes_key}"] = {
            'qq': qq_fig.to_dict()
        }

    # convert plot data to JSON-like structure
    plot_data_json = json.dumps(plot_data)

    #
    # Global dNdS Plot

    mle = df_global.mle
    upper_global = df_global.cihigh - mle
    lower_global = mle - df_global.cilow

    dict_erry_global = dict(
        type='data',
        symmetric=False,
        array=np.round(upper_global, 2).tolist(),
        arrayminus=np.round(lower_global, 2).tolist(),
        thickness=thk_err_global,
        width=wid_err_global,
        color=col_err_nonsig
    )

    x_labels = ['Missense', 'Nonsense', 'Splice site', 'Truncating', 'All']
    x_global = [i + 1 for i in range(len(df_global))]
    xlim_global = [-hor_buffer_global + min(x_global), max(x_global) + hor_buffer_global]

    fig_dnds_global = go.Figure()
    fig_dnds_global.add_trace(
        go.Scatter(
            x=xlim_global,
            y=[1, 1],
            mode='lines',
            line=dict(dash=typ_thick, color=col_thick, width=thk_thick),
            showlegend=False,
            hoverinfo='skip'
        )
    )
    fig_dnds_global.add_trace(
        go.Scatter(
            x=x_global,
            y=df_global['mle'].tolist(),
            error_y=dict_erry_global,
            mode='markers',
            marker=dict(color=col_nonsig, opacity=opac_nonsig),
            showlegend=False,
            xhoverformat='.2f',
            yhoverformat='.2f'
        )
    )

    fig_dnds_global.update_layout(
        title=dict(
            text='dNdS' + ' ratio across all mutations:'.title(),
            # font=dict(size=20)  # Set your desired title font size
        ),
        xaxis_title='',
        yaxis_title='dNdS ratio',
        xaxis=dict(
            range=xlim_global,
            tickvals=x_global,
            ticktext=x_labels,
            tickfont=dict(size=fs_axes)  # size of the tick labels
        ),
        yaxis=dict(
            titlefont=dict(size=fs_axes)  # size of y-axis title
        ),
        template='plotly_white'
    )

    #
    # Missense dNdS Plot

    fig_dnds_mis = plot_dnds(df_merged, 'mis', 'missense', alp, alp_nearsig)

    #
    # Truncating dNdS Plot

    fig_dnds_tru = plot_dnds(df_merged, 'tru', 'truncating', alp, alp_nearsig)

    #
    # Table Plot

    n_rows = max(n_rows_min, int(np.sum(ind_sig) * (1 + n_rows_buffer)))
    df_plot = df_plot.iloc[:n_rows][cols_kept].copy()

    for col in ['PVAL', 'FDR'] + list(dict_pvals.values()):
    # for col in ['PVAL', 'FDR']:
        df_plot[col] = reformat_numbers(df_plot[col].to_numpy())
    for col in list(dict_w.values()):
        df_plot[col] = reformat_numbers(df_plot[col].to_numpy(), format='{:.2f}')
    for c in list(dict_n.values()) + ['RANK']:
        df_plot[c] = df_plot[c].astype(int)

    headerColor = 'grey'
    rowEvenColor = 'lightgrey'
    rowOddColor = 'white'
    # making the significant rows bold
    df_plot = df_plot.astype(str)
    for i in range(df_plot.shape[0]):
        if ind_sig[i]:
            df_plot.loc[i, :] = '<b>' + df_plot.loc[i, :].astype(str) + '</b>'
    # generate table figure
    fig_table = go.Figure(data=[go.Table(
        header=dict(values=['<b>' + col + '</b>' for col in df_plot.columns],
                    line_color='darkslategray',
                    fill_color=headerColor,
                    align=['left'] + ['center'] * (len(df_plot.columns) - 1),
                    font=dict(color='white', size=12)
                    ),
        cells=dict(values=[df_plot[col].tolist() for col in df_plot.columns],
                   line_color='darkslategray',
                   fill_color=[[rowOddColor if i % 2 == 0 else rowEvenColor for i in range(df_plot.shape[0])]],
                   align=['left'] + ['center'] * (len(df_plot.columns) - 1),
                   font=dict(color='darkslategray', size=11),
                   # format=['html'] * len(df_plot.columns)  # Enable HTML formatting
                   )
    )
    ])

    #
    # Assemble Final HTML

    # save static figures as html divs
    fig_dnds_global_html = fig_dnds_global.to_html(full_html=False, include_plotlyjs='cdn')
    fig_dnds_mis_html = fig_dnds_mis.to_html(full_html=False, include_plotlyjs='cdn')
    fig_dnds_tru_html = fig_dnds_tru.to_html(full_html=False, include_plotlyjs='cdn')
    fig_table_html = fig_table.to_html(full_html=False, include_plotlyjs='cdn')

    # combine everything into the final html
    html_content = html_content.format(
        ngenes_options=ngenes_options,
        plot_data=plot_data_json,
        fig_dnds_global_html = fig_dnds_global_html,
        fig_dnds_mis_html=fig_dnds_mis_html,
        fig_dnds_tru_html=fig_dnds_tru_html,
        fig_table_html=fig_table_html
    )

    # save to an html file
    with open(dir_output + '/' + ('' if (prefix_output is None) else prefix_output + '_') +
              'dndscv_report.html', 'w') as f:
        f.write(html_content)


def gen_dndscv_results_app_component(
    dnd_df,
    num_gene_values,
    dnd_radio_item_selection

):
    """ 
    """
    debugging= ""
    all_page_content = []
    # dig_df = data.df[SNV_DATA_COLUMN_NAME][0][DIG_DATAFRAME_IDX].copy() # gets the dig report data for a specific cohort
    dnd_df = dnd_df.copy()
    figure1 = go.Figure()
    all_page_content = [
            figure1,
            dnd_radio_item_selection,
            RESULTS_DROPDOWN_VALUES
        ]

    # NEED TO ASK DAVID WHICH COLUMN TO SORT THE dnd_df by
    # dnd_df = dnd_df.sort_values(by='PVAL'+ "_"+ combined_mutation_type[mutation_type] + "_" + scatterpoint_type[p_val_type])
    # dnd_df['RANK'] = np.array([i+1 for i in range(len(dnd_df))])

    return all_page_content