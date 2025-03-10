import plotly.graph_objects as go
import pandas as pd
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import numpy as np
from SNVReviewers.AppComponents.utils import reformat_numbers

MUTSIG_REPORT_COLUMN_NAMES = ["rank", "gene", "nnei", "nncd", "nsil", "nmis", "nstp", "nspl", "nind",
                              "nnon", "npat", "pCV", "pCL", "pFN", "pCL2", "pFN2", "p", 
                            ]
                            # ["RANK", "GENE", "NNEI", "NNCD", "NSIL", "NMIS", "NSTP", "NSPL", "NIND",
                            #   "NNON", "NPAT", "PCV", "PCL", "PFN", "PCL2", "PFN2", "P", 
                            # #   "FDR", "CGC", "PANCAN"
                            # ]

def gen_mutsig_results_app_component(
    df_mutsig,
    num_gene_values,
):
    """ 
    """
    debugging= ""
    all_page_content = []
    df_mutsig = df_mutsig.copy()
    results_table = go.Figure()
    qq_fig = go.Figure()

    # df_mutsig = df_mutsig.sort_values(by='fdr')
    df_mutsig['rank'] = np.array([i+1 for i in range(len(df_mutsig))])

    # qq_fig, df_plot = generate_mutsig_report(mutsig_df_plot, num_gene_values)

    for clm_nm in MUTSIG_REPORT_COLUMN_NAMES:
        
        # reformats the columns with float values in them
        if pd.api.types.is_float_dtype(df_mutsig[clm_nm]) and clm_nm != 'rank':
            df_mutsig[clm_nm] = reformat_numbers(df_mutsig[clm_nm], format='{:.3E}')
    
    # NEED TO REFORMAT THE dnd_df_plot dataframe for 2 point decimal precision
    all_page_content = [
            # data to be displayed in the tables for the results dndscv report
            df_mutsig.to_dict('records'),

            # result mutsig report figures
            qq_fig,
            
            # mutsig dropdown menu value
            num_gene_values,

            debugging
        ]

    return all_page_content