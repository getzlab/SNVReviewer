def gen_dig_combined_app_component_data_internal_callback(
    data: GenericData,
    idx,
    dig_label,
    dig_type_selection, # radio item selection
    mutation_type,
    burden_type,
    p_val_type,
    display_toggle_value,
    display_label_value
):
    """
    
    """
    all_page_content = []

    # check if mutation_type is an old dropdown (coding region, non coding regions) menu option
    if mutation_type == "":
        mutation_type = 'indels_snvs'
        burden_type = 'total'
        p_val_type = 'uniform_p_mid'

    dig_df = data.df[SNV_DATA_COLUMN_NAME][0][DIG_DATAFRAME_IDX].copy() # gets the dig report data for a specific cohort

    # if dig_label == "Combined":
    dig_df = dig_df.sort_values(by='PVAL'+ "_"+ combined_mutation_type[mutation_type] + "_" + scatterpoint_type[p_val_type])
    dig_df['RANK'] = np.array([i+1 for i in range(len(dig_df))])

    debugging_component = ""
    dig_output_type = "Combined"
    display_bounds = display_toggle_value # whether to display the bounds
    
    # dir_output = "./example_notebooks/data/dig_data"
    qq_fig = go.Figure()
    fig_mu = go.Figure()
    fig_sigma = go.Figure()
    dnds_fig = go.Figure()
    volcano_fig = go.Figure()

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
        
    # MAKE SURE TO REDO THE CODING REGION VALUES, ADD IF STATEMENTS TO CHANGE WHAT GETS DISPLAYED BASED ON THE DROP DOWN MENU
    # volcano_fig, qq_fig, fig_mu, fig_sigma, dnds_fig, table_fig, text_special = generate_coding_region_report(dig_df, mutation_type, burden_type, display_labels_key, p_val_type)

    dig_data_columns = []
    
    for clm_nm in DIG_REPORT_COMBINED_COLUMN_NAMES:
        dig_data_clm_dict = {}
        dig_data_clm_dict["name"] = clm_nm
        dig_data_clm_dict["id"] = clm_nm

        if clm_nm == 'FDR' or "PVAL" in clm_nm:
            # gets the column name corresponding to the mutation type and the scatterpoint type
            dig_data_clm_dict["id"] = clm_nm + "_"+ combined_mutation_type[mutation_type] + "_" + scatterpoint_type[p_val_type]

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
        
    # ONLY GETTING THE FIRST 100 ROWS OF DATA TO DISPLAY IN THE TABLE
    # REMOVE THE DEBUGGING LATER!!!   
    dig_df = dig_df[:100]

    all_page_content = [
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
                COMBINED_MUT_DROPDOWN,
                COMBINED_BUR_DROPDOWN,
                COMBINED_SCATTER_DROPDOWN, 
                text_special,
                debugging_component
            ]
        
    return all_page_content