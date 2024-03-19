#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import networkx as nx

"""
Published on Sat Mar 2 19:43:26 2024

@author: Eric Li
"""
"""
This code obtains closeness centrality scores of all areas via connection matrices and
correlation plots with hierarchy scores.
"""


# # Set input/output directories

# In[ ]:


input_dir = r'./Input/'
output_dir = r'./Output/'
results_dir = r'./results/'


# # Obtain hierarchy results

# In[ ]:


xls=pd.ExcelFile(results_dir+"hierarchy_summary_antero_retro_CreConf.xlsx")
df = pd.read_excel(xls, "hierarchy_all_regions")

hierarchy = pd.concat([df["areas"], df["CC+TCCT iterated"]], axis = 1)
hierarchy = hierarchy.values


# # Obtain + preprocess normalized connection matrices. Set NaN values to zero and transpose retrograde/combined matrices.

# In[ ]:


xls1=pd.ExcelFile(input_dir+"combined_ipsi_VN.xlsx")

df1=pd.read_excel(xls1,'combined_ipsi_VN')
df1 = df1.set_index("area")

xls2=pd.ExcelFile(input_dir+"ipsi_C57_antero_VN.xlsx")

df2=pd.read_excel(xls2,'ipsi_C57_antero_VN')
df2 = df2.iloc[:, 1:]
df2 = df2.set_index("source_area")

xls3=pd.ExcelFile(input_dir+"retro_ipsi_VN.xlsx")

df3=pd.read_excel(xls3,'retro_ipsi_VN')
df3 = df3.iloc[:, 1:]
df3 = df3.set_index("target")

areas = df1.columns.values

matrix_nan = np.isnan(df1.values)
antero_nan = np.isnan(df2.values)
retro_nan = np.isnan(df3.values)

df1.iloc[matrix_nan] = 0
df2.iloc[antero_nan] = 0
df3.iloc[retro_nan] = 0

matrix = df1.values
matrix = np.transpose(matrix)
antero = df2.values
retro = df3.values
retro = np.transpose(retro)


# # Define functions needed

# In[ ]:


def get_module(area):
    if (area == 'GU' or area == 'VISC' or area == 'SSs' or 
        area == 'SSp-bfd' or area == 'SSp-ll' or 
        area == 'SSp-ul' or area == 'SSp-un' or 
        area == 'SSp-n' or area == 'SSp-m' or 
        area == 'MOp' or area == 'MOs'):
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
             area == 'SSp-tr' or area == 'RSPagl' or 
             area == 'RSPd' or area == 'RSPv'):
        return 'blue'
    elif (area == 'TEa' or area == 'AUDd' or 
             area == 'AUDp' or area == 'AUDpo' or 
             area == 'AUDv'):
        return 'purple'
    else:
        return 'black'

def make_plot(centrality, name):
    hierarchy_values = []
    centrality_values = []
    centrality_areas = []
    f = plt.figure()  
    for key, item in centrality.items():
        for i in range(len(hierarchy)):
            if areas[key] == hierarchy[i, 0]:
                centrality_areas.append(areas[key])
                centrality_values.append(item)
                hierarchy_values.append(hierarchy[i, 1])
                plt.text(item, hierarchy[i, 1], hierarchy[i, 0])
                plt.plot(item, hierarchy[i, 1],'o', color = get_module(hierarchy[i, 0]))
                break
          
    plt.xlabel('centrality value')
    plt.ylabel('hierarchy score')
    plt.title(name+', r='+str(np.corrcoef(centrality_values, hierarchy_values)[0, 1]))
    plt.show()
    #f.savefig(output_dir+name+".svg", format="svg")
    return centrality_areas, centrality_values, hierarchy_values

def create_networkx(matrix, weight, in_out):
    G = nx.DiGraph()

    num_nodes = matrix.shape[0]
    G.add_nodes_from(range(num_nodes))

    for i in range(num_nodes):
        for j in range(num_nodes):
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

# In[ ]:


antero_4 = sparsify(antero, None, 10**(-4), 'v')
antero_binary = create_networkx(antero_4, False, 'in')
centrality_areas, antero_binary_centrality, antero_hierarchy_scores = make_plot(antero_binary, 'anterograde, inward centrality, binarized, threshold @ 10^-4')

antero_data = [centrality_areas, antero_binary_centrality, antero_hierarchy_scores]
antero_data = np.transpose(antero_data)
columns = ['areas', 'anterograde (inward) centrality', 'hierarchy scores']
df = pd.DataFrame(antero_data, columns=columns)
df = df.sort_values('anterograde (inward) centrality', ascending = False)
#df.to_excel(output_dir+'antero_inward_centrality.xlsx')

retro_4 = sparsify(retro, None, 10**(-4), 'v')
retro_binary = create_networkx(retro_4, False, 'out')
centrality_areas, retro_binary_centrality, retro_hierarchy_scores = make_plot(retro_binary, 'retrograde, outward centrality, binarized, threshold @ 10^-4')

retro_data = [centrality_areas, retro_binary_centrality, retro_hierarchy_scores]
retro_data = np.transpose(retro_data)
columns = ['areas', 'retrograde (outward) centrality', 'hierarchy scores']
df = pd.DataFrame(retro_data, columns=columns)
df = df.sort_values('retrograde (outward) centrality', ascending = False)
#df.to_excel(output_dir+'retro_outward_centrality.xlsx')

averaged = combine_dictionaries(antero_binary, retro_binary)
centrality_areas, averaged_centrality, averaged_hierarchy_scores = make_plot(averaged, 'averaged, antero (inward) + retro (outward) centrality, binarized, threshold @ 10^-4')

averaged_data = [centrality_areas, averaged_centrality, averaged_hierarchy_scores]
averaged_data = np.transpose(averaged_data)
columns = ['areas', 'antero (inward) + retro (outward) centrality', 'hierarchy scores']
df = pd.DataFrame(averaged_data, columns=columns)
df = df.sort_values('antero (inward) + retro (outward) centrality', ascending = False)
#df.to_excel(output_dir+'averaged_centrality.xlsx')

