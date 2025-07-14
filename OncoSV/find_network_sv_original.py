#!/usr/bin/env python3

import pandas as pd
import networkx as nx

def group_by_group(df):
    # Sort the DataFrame alphabetically by the 'ID' column
    df = df.sort_values(by='ID')

    # Splitting and getting unique values from 'ID' column
    sv_list = set(df['ID'].str.split(',').sum())

    complex_list = [{'ID': None, 'group': None, 'CHROM': None, 'CHROM2': None, 'POS': None, 'END': None, 'Read_Count': None, 'SV_Count': None, 'AF': None, 'CSV_Type': None, 'POS_BKPT': None, 'END_BKPT': None, 'final_combination': None, 'any_overlapping_sv': None}]
    group_number = 0

    for i in sv_list:
        # Use regex to ensure exact match of i within the 'ID' strings
        row_i = df.index[df['ID'].str.contains(f"\\b{i}\\b", regex=True)].tolist()

        if len(row_i) > 1:
            sv_i = df.loc[row_i, 'ID'].tolist()
            edges = []
            detected = [False] * len(sv_i)

            # Convert each ID string into a list of identifiers to maintain order
            sv_lists = [item.split(',') for item in sv_i]

            for j in range(len(sv_i)):
                for z in range(j + 1, len(sv_i)):  # Avoid repeating pairs and self-comparison
                    if sv_lists[j] == sv_lists[z]:
                        edges.append((j, z))
                        detected[j] = True
                        detected[z] = True

            for k in range(len(sv_i)):
                if not detected[k]:
                    edges.append((k, k))

            if edges:
                group_number += 1
                G = nx.Graph()

                # Add valid edges without weights
                for edge in edges:
                    if len(edge) == 2:
                        G.add_edge(*edge)

                groups = nx.connected_components(G)
                grouped_nodes = {idx: list(comp) for idx, comp in enumerate(groups)}

                for node_i, grouped_i in grouped_nodes.items():
                    if len(grouped_i) > 1:
                        group_row = [sv_i[idx] for idx in grouped_i]
                        number = max(range(len(group_row)), key=lambda idx: group_row[idx].count(','))
                        selected = group_row[number]
                    else:
                        selected = sv_i[grouped_i[0]]

                    existing_index = [idx for idx, x in enumerate(complex_list) if x['ID'] == selected]

                    if len(existing_index) == 1:
                        complex_list[existing_index[0]]['group'] = f"{complex_list[existing_index[0]]['group']};{group_number}"
                    else:
                        # Find the corresponding values from the DataFrame
                        row_data = df.loc[df['ID'] == selected].iloc[0]
                        complex_list.append({
                            'ID': selected,
                            'group': group_number,
                            'CHROM': row_data['CHROM'],
                            'CHROM2': row_data['CHROM2'],
                            'POS': row_data['POS'],
                            'END': row_data['END'],
                            'Read_Count': row_data['Read_Count'],
                            'SV_Count': row_data['SV_Count'],
                            'AF': row_data['AF'],
                            'CSV_Type': row_data['CSV_Type'],
                            'POS_BKPT': row_data['POS_BKPT'],
                            'END_BKPT': row_data['END_BKPT'],
                            'final_combination': row_data['final_combination'],
                            'any_overlapping_sv': row_data['any_overlapping_sv']
                        })
        elif len(row_i) == 1:
            selected = df.loc[row_i[0], 'ID']
            existing_index = [idx for idx, x in enumerate(complex_list) if x['ID'] == selected]

            if len(existing_index) == 1:
                continue

            group_number += 1
            # Find the corresponding values from the DataFrame
            row_data = df.loc[df['ID'] == selected].iloc[0]
            complex_list.append({
                'ID': selected,
                'group': group_number,
                'CHROM': row_data['CHROM'],
                'CHROM2': row_data['CHROM2'],
                'POS': row_data['POS'],
                'END': row_data['END'],
                'Read_Count': row_data['Read_Count'],
                'SV_Count': row_data['SV_Count'],
                'AF': row_data['AF'],
                'CSV_Type': row_data['CSV_Type'],
                'POS_BKPT': row_data['POS_BKPT'],
                'END_BKPT': row_data['END_BKPT'],
                'final_combination': row_data['final_combination'],
                'any_overlapping_sv': row_data['any_overlapping_sv']
            })
            df.loc[row_i, 'group'] = group_number

    complex_list = [x for x in complex_list if x['ID'] is not None]
    complex_df = pd.DataFrame(complex_list)

    return complex_df

class UnionFind:
    def __init__(self):
        self.parent = {}

    def find(self, item):
        if item not in self.parent:
            self.parent[item] = item
            return item
        if self.parent[item] == item:
            return item
        self.parent[item] = self.find(self.parent[item])  # Path compression
        return self.parent[item]

    def union(self, set1, set2):
        root1 = self.find(set1)
        root2 = self.find(set2)
        if root1 != root2:
            self.parent[root2] = root1  # Union

def identify_networks(df):
    uf = UnionFind()
    # Parsing and union of group IDs
    for index, row in df.iterrows():
        groups = str(row['group']).split(';')
        base = groups[0]
        for m in groups:
            uf.union(base, m)

    # Assign Network IDs
    network_ids = {m: idx + 1 for idx, m in enumerate(set(uf.find(m) for m in uf.parent))}

    # Add Network IDs to the DataFrame
    df['Network'] = df['group'].apply(lambda x: network_ids[uf.find(str(x).split(';')[0])])

    return df

