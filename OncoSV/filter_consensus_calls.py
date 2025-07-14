#!/usr/bin/env python3

import pandas as pd
import re

def filter_consensus_calls(df):
    # Extract 'sv_caller' values from the 'ID' column
    df['ID'] = df['ID'].astype(str)
    df['sv_caller'] = [re.search(r"^(.*?)\.", item).group(1) for item in df['ID']]

    # Calculate the number of unique 'sv_caller' values for each 'ConsensusSV_ID'
    NUM_CALLERS = df.groupby('ConsensusSV_ID')['sv_caller'].nunique()

    # Add 'consensus_tool_count' to the DataFrame
    df['NUM_CALLERS'] = df['ConsensusSV_ID'].map(NUM_CALLERS)

    # Create a new DataFrame with duplicated 'ConsensusSV_ID' values
    duplicated_consensus_ids = NUM_CALLERS[NUM_CALLERS >= 2].index
    df_select = df[df['ConsensusSV_ID'].isin(duplicated_consensus_ids)]

    # Create a function to prioritize Sniffles2 and cuteSV and return the index
    def prioritise_tools(group):
        if 'Sniffles2' in group['sv_caller'].values:
            return group[group['sv_caller'] == 'Sniffles2'].index
        else:
            return group[group['sv_caller'] == 'cuteSV'].index

    # Apply the prioritize_tools function and create a Series with the same index
    priority_series = df_select.groupby('ConsensusSV_ID').apply(prioritise_tools)
    index_numbers_list = [index[0] for index in priority_series]

    # Extract the rows with the priority tool based on the index
    filtered_df = df_select.loc[index_numbers_list]
    # Remove unnecessary columns
    filtered_df = filtered_df.drop(columns=['sv_caller'])

    #sort the table based on CHROM and process_sample_data
    filtered_df = filtered_df.sort_values(by=['CHROM', 'POS'])

    # Replace "type" with the values from the SVTYPE column
    filtered_df['ConsensusSV_ID'] = filtered_df.apply(lambda row: \
    row['ConsensusSV_ID'].replace('type', row['SVTYPE']), axis=1)

    return filtered_df
