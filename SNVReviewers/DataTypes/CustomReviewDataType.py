from AnnoMate.DataTypes.GenericData import GenericData
from AnnoMate.Data import Data, DataAnnotation

import pandas as pd
import numpy as np
from typing import Union, List, Dict
from pathlib import Path
import os

class CustomReviewDataType(Data):
    def __init__(
        self,
        index,
        description,
        df: pd.DataFrame,
        annot_df: pd.DataFrame = None,
        annot_col_config_dict: Dict = None,
        history_df: pd.DataFrame = None,
    ):
        
        super().__init__(
            index=index,
            description=description,
            annot_df=annot_df,
            annot_col_config_dict=annot_col_config_dict,
            history_df=history_df,
        )
        
        # custom formatting of df
        
        self.df = df