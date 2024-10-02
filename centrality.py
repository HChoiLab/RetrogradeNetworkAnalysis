#!/usr/bin/env python
# coding: utf-8


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import networkx as nx

"""
Published on Tues Oct 1 19:43:26 2024

@author: Eric Li
"""
"""
This code obtains closeness centrality scores of all areas via connection matrices and
correlation plots with hierarchy scores.
"""


# # Set input/output directories


input_dir = r'./Input/'
output_dir = r'./Output/'
results_dir = r'./results/'


# # Obtain hierarchy results


xls=pd.ExcelFile(results_dir+"hierarchy_summary.xlsx")
df = pd.read_excel(xls, "hierarchy_all_regions")

hierarchy = pd.concat([df["areas"], df["CC+TCCT+CB iterated"]], axis = 1)
hierarchy = hierarchy.values
hierarchy = np.concatenate([hierarchy[:14],hierarchy[22:36], hierarchy[37:]], axis=0) #remove CP areas and OT_9


# # Obtain + preprocess normalized connection matrices. Set NaN values to zero and transpose retrograde matrix.

def reorder_matrix(df):
    m = df.shape[0]
    n = df.shape[1]
    
    row_names = df.index.tolist()
    column_names = df.columns.tolist()
    
    sorted_row_names = sorted(row_names)
    sorted_column_names = sorted(column_names)
    
    if m < n:
        remaining_columns = [col for col in sorted_column_names if col not in sorted_row_names]
        sorted_column_names = sorted_row_names + remaining_columns
        areas = sorted_column_names
    elif m > n:
        remaining_rows = [row for row in sorted_row_names if row not in sorted_col_names]
        sorted_row_names = sorted_column_names + remaining_rows
        areas = sorted_row_names
    else:
        areas = sorted_row_names
        
    df = df.reindex(index=sorted_row_names, columns=sorted_column_names)
    return df, areas
    
xls1=pd.ExcelFile(input_dir+"antero_ipsi_VN.xlsx")

df1=pd.read_excel(xls1,'antero_ipsi_VN')
df1 = df1.iloc[:, 1:]
df1 = df1.set_index("source_area")

xls2=pd.ExcelFile(input_dir+"retro_ipsi_VN.xlsx")

df2=pd.read_excel(xls2,'retro_ipsi_VN')
df2 = df2.iloc[:, 1:]
df2 = df2.set_index("target")

antero_nan = np.isnan(df1.values)
retro_nan = np.isnan(df2.values)

if len(df1.iloc[antero_nan]) != 0:
    df1.iloc[antero_nan] = 0
if len(df2.iloc[retro_nan]) != 0:
    df2.iloc[retro_nan] = 0

df1, areas1 = reorder_matrix(df1)
df2, areas2 = reorder_matrix(df2)

antero = df1.values
retro = df2.values
retro = np.transpose(retro)


# # Define functions needed


def get_module(area):
    if (area == 'GU' or area == 'VISC' or area == 'SSs' or 
        area == 'SSp-bfd' or area == 'SSp-tr' or
        area == 'SSp-ll' or area == 'SSp-ul' or
        area == 'SSp-un' or area == 'SSp-n' or
        area == 'SSp-m' or area == 'MOp' or
        area == 'MOs'):
        return 'orange'
    elif (area == 'FRP' or area == 'PL' or 
             area == 'ILA' or area == 'ORBl' or 
             area == 'ORBm' or area == 'ORBvl' or 
             area == 'AId' or area == 'AIv' or 
             area == 'AIp' or area == 'ECT'):
        return 'red'
    elif (area == 'VISal' or area == 'VISl' or 
             area == 'VISp' or area == 'VISpl' or 
             area == 'VISli' or area == 'VISpor' or 
             area == 'VISrl' or area == 'VISa' or 
             area == 'VISam' or area == 'VISpm'):
        return 'cyan'
    elif (area == 'ACAd' or area == 'ACAv' or 
          area == 'RSPagl' or area == 'RSPd' or 
          area == 'RSPv'):
        return 'blue'
    elif (area == 'TEa' or area == 'AUDd' or 
             area == 'AUDp' or area == 'AUDpo' or 
             area == 'AUDv'):
        return 'purple'
    else:
        return 'black'

def make_plot(centrality, hierarchy, name):
    f = plt.figure()  
    for area, hierarchy_value in hierarchy.items():
        plt.text(centrality[area], hierarchy_value, area)
        plt.plot(centrality[area], hierarchy_value,'o', color = get_module(area))

    areas = np.array(list(hierarchy.keys()))
    centrality_array = np.array(list(centrality.values())).astype(float)
    hierarchy_array = np.array(list(hierarchy.values())).astype(float)
    plt.xlabel('centrality value')
    plt.ylabel('hierarchy score')
    plt.title(name+', r='+str(np.corrcoef(centrality_array, hierarchy_array)[0, 1]))
    plt.show()

    return np.transpose([areas, centrality_array, hierarchy_array])

def get_hierarchy_areas(areas, centrality, hierarchy):
    centrality_areas = []
    centrality_scores = []
    hierarchy_scores = []
    for i in range(len(hierarchy)):
        for j in range(len(areas)):
            if areas[j] == hierarchy[i, 0]:
                centrality_areas.append(areas[j])
                centrality_scores.append(centrality[j])
                hierarchy_scores.append(hierarchy[i, 1])
    return dict(zip(centrality_areas, centrality_scores)), dict(zip(centrality_areas, hierarchy_scores))

def create_networkx(matrix, weight, in_out):
    G = nx.DiGraph()

    num_sources = matrix.shape[0]
    num_targets = matrix.shape[1]
    if num_sources > num_targets:
        G.add_nodes_from(range(num_sources))
    elif num_sources < num_targets:
        G.add_nodes_from(range(num_targets))
    else:
        G.add_nodes_from(range(num_sources))

    for i in range(num_sources):
        for j in range(num_targets):
            if matrix[i][j] != 0:
                G.add_edge(i, j, weight=matrix[i][j])
    if in_out == 'out':
        G = G.reverse()
    if weight:
        centrality = nx.closeness_centrality(G, distance="weight", wf_improved=True)
    else:
        centrality = nx.closeness_centrality(G, wf_improved=True)
    return centrality

def sparsify(matrix, percentage, value, pv):
    matrix_sorted = np.sort(matrix.flatten())
    if pv == 'p':
        threshold_index = int((1 - percentage) * len(matrix_sorted))
        matrix_threshold = matrix_sorted[threshold_index]
    else:
        matrix_threshold = value
        
    matrix_binary = np.where(matrix > matrix_threshold, 1, 0)
    return matrix_binary

def combine_dictionaries(dict1, dict2):
    combined_dict = {}

    for key in dict1.keys():
        if key in dict2:
            average_value = (dict1[key] + dict2[key]) / 2
            combined_dict[key] = average_value
        else:
            combined_dict[key] = dict1[key]

    for key in dict2.keys():
        if key not in dict1:
            combined_dict[key] = dict2[key]

    return combined_dict

# # Obtain closeness centrality values and make correlation plots with hierarchy scores.


antero_4 = sparsify(antero, None, 10**(-4), 'v')
antero_centrality_all = create_networkx(antero_4, False, 'in')
antero_centrality, antero_hierarchy = get_hierarchy_areas(areas1, antero_centrality_all, hierarchy)
antero_data = make_plot(antero_centrality, antero_hierarchy, 'anterograde, inward centrality, binarized, threshold @ 10^-4')

columns = ['areas', 'anterograde (inward) centrality', 'hierarchy scores']
df = pd.DataFrame(antero_data, columns=columns)
df = df.sort_values('anterograde (inward) centrality', ascending = False)
df.to_excel(output_dir+'antero_inward_centrality.xlsx')

retro_4 = sparsify(retro, None, 10**(-4), 'v')
retro_centrality_all = create_networkx(retro_4, False, 'out')
retro_centrality, retro_hierarchy = get_hierarchy_areas(areas2, retro_centrality_all, hierarchy)
retro_data = make_plot(retro_centrality, retro_hierarchy, 'retrograde, outward centrality, binarized, threshold @ 10^-4')

columns = ['areas', 'retrograde (outward) centrality', 'hierarchy scores']
df = pd.DataFrame(retro_data, columns=columns)
df = df.sort_values('retrograde (outward) centrality', ascending = False)
df.to_excel(output_dir+'retro_outward_centrality.xlsx')

averaged_centrality = combine_dictionaries(antero_centrality, retro_centrality)
averaged_hierarchy = combine_dictionaries(antero_hierarchy, retro_hierarchy)
averaged_data = make_plot(averaged_centrality, averaged_hierarchy, 
                          'averaged, antero (inward) + retro (outward) centrality, binarized, threshold @ 10^-4')

columns = ['areas', 'antero (inward) + retro (outward) centrality', 'hierarchy scores']
df = pd.DataFrame(averaged_data, columns=columns)
df = df.sort_values('antero (inward) + retro (outward) centrality', ascending = False)
df.to_excel(output_dir+'averaged_centrality.xlsx')