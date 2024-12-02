import pandas as pd
from typing import List, Dict

from dash.dependencies import State
from rpy2.robjects import pandas2ri

from MyCustomReviewers.AppComponents.ACustomAppComponent import gen_custom_app_component

from AnnoMate.Data import DataAnnotation
from AnnoMate.ReviewDataApp import ReviewDataApp
from AnnoMate.DataTypes.GenericData import GenericData
from AnnoMate.ReviewerTemplate import ReviewerTemplate
from AnnoMate.AppComponents.DataTableComponents import gen_annotated_data_info_table_component
from AnnoMate.AnnotationDisplayComponent import NumberAnnotationDisplay


class SNVReviewer(ReviewerTemplate):
    def gen_data(self,
                 description: str,
                 df: pd.DataFrame,
                 index: List,
                 annot_df: pd.DataFrame = None,
                 annot_col_config_dict: Dict = None,
                 history_df: pd.DataFrame = None) -> GenericData:
        """Set up data for CustomReviewer
        Returns
        -------
        GenericData
            Data object that contains only one dataframe
        """

        return GenericData(index=index,
                           description=description,
                           df=df,
                           annot_df=annot_df,
                           annot_col_config_dict=annot_col_config_dict,
                           history_df=history_df)

    def gen_review_app(self,
                       sample_info_cols,
                       acs_col,
                       csize=None,
                       step_size=None
                       ) -> ReviewDataApp:
        """
        Parameters
        ==========
        sample_info_cols: list[str]
            List of columns in data
        acs_col: str
            Column name in data with path to seg file from alleliccapseg
        csize: dict
            Dictionary with chromosome sizes
        step_size: float
            Minimum increment allowed for purity (default is 0.01)
        """
        app = ReviewDataApp()
        
        app.add_component(
            gen_custom_app_component(), 
        )

        return app

    def set_default_autofill(self):
        """
        self.add_autofill(<button name>, State('<id in layout>', 'attribute of object'), <annot name>)
        """
        pass

    def set_default_review_data_annotations(self):
        """
        self.add_review_data_annotation(
            annot_name='annot name',
            review_data_annotation=DataAnnotation('float')
        )
        """
        pass


    def set_default_review_data_annotations_app_display(self):
        """
        self.add_annotation_display_component(<annot name>, NumberAnnotationDisplay())
        """
        pass
        
        