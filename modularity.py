#!/usr/bin/env python
# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import bct as bct

"""
Published on Sat Mar 2 19:03:31 2024

@author: Eric Li
"""
"""
This code generates modules of the anterograde, retrograde, and averaged cortical connection 
matrices and compares their modularity scores to those of shuffled versions.
"""


# # Set input and output directories


input_dir = r'./Input/'
output_dir = r'./Output/'


# # Read and preprocess excel file with connection matrices. Set all NaN values to zero and obtain cortical areas. Transpose retrograde and combined matrices.


xls1=pd.ExcelFile(input_dir+"combined_ipsi.xlsx")

df1=pd.read_excel(xls1,'in')
df1 = df1.set_index("target")

xls2=pd.ExcelFile(input_dir+"antero_ipsi.xlsx")

df2=pd.read_excel(xls2,'antero_ipsi')
df2 = df2.iloc[:, 1:]
df2 = df2.set_index("source_area")

xls3=pd.ExcelFile(input_dir+"retro_ipsi.xlsx")

df3=pd.read_excel(xls3,'retro_ipsi')
df3 = df3.set_index("target")

df1 = df1.drop("PERI", axis = 0)
df1 = df1.drop("PERI", axis = 1)
df2 = df2.drop("PERI", axis = 0)
df2 = df2.drop("PERI", axis = 1)
df3 = df3.drop("PERI", axis = 0)
df3 = df3.drop("PERI", axis = 1)

areas = df1.columns.values

cortical_areas = areas[:42]
matrix_cortical = df1.iloc[:42, :42]
antero_cortical = df2.iloc[:42, :42]
retro_cortical = df3.iloc[:42, :42]

matrix_nan = np.isnan(matrix_cortical.values)
antero_nan = np.isnan(antero_cortical.values)
retro_nan = np.isnan(retro_cortical.values)

matrix_cortical.iloc[matrix_nan] = 0
antero_cortical.iloc[antero_nan] = 0
retro_cortical.iloc[retro_nan] = 0

m = matrix_cortical.values
m = np.transpose(m)
a = antero_cortical.values
r = retro_cortical.values
r = np.transpose(r)


# # Create shuffled version of cortical connection matrices


matrix_cortical_shuffled = matrix_cortical.sample(frac=1).reset_index(drop=True)
matrix_cortical_shuffled = matrix_cortical_shuffled.sample(frac=1, axis=1).reset_index(drop=True)

antero_cortical_shuffled = antero_cortical.sample(frac=1).reset_index(drop=True)
antero_cortical_shuffled = antero_cortical_shuffled.sample(frac=1, axis=1).reset_index(drop=True)

retro_cortical_shuffled = retro_cortical.sample(frac=1).reset_index(drop=True)
retro_cortical_shuffled = retro_cortical_shuffled.sample(frac=1, axis=1).reset_index(drop=True)

m_shuffled = matrix_cortical_shuffled.values
m_shuffled = np.transpose(m_shuffled)
a_shuffled = antero_cortical_shuffled.values
r_shuffled = retro_cortical_shuffled.values
r_shuffled = np.transpose(r_shuffled)


# # Obtain communities using Louvain Algorithm and compare to shuffled matrices.


name = 'retrograde modularity'

num = 16
data_mod_scores = np.zeros([num, 2])
shuffled_mod_scores = np.zeros([num, 2])
g = 0
for i in range(num):
    [[community], mod_score] = bct.modularity_louvain_dir(r, gamma=g, hierarchy=True, seed=1234)
    [[community_shuffled], mod_score_shuffled] = bct.modularity_louvain_dir(r_shuffled, gamma=g, hierarchy=True)
    data_mod_scores[i, 0] = i*0.1
    data_mod_scores[i, 1] = mod_score[0]
    shuffled_mod_scores[i, 0] = i*0.1
    shuffled_mod_scores[i, 1] = mod_score_shuffled[0]
    g += 0.1
    
max_gamma = np.argmax(data_mod_scores[:, 1] - shuffled_mod_scores[:, 1])
max_gamma_value = np.max(data_mod_scores[:, 1] - shuffled_mod_scores[:, 1])   

fig, ax = plt.subplots()
ax.plot(data_mod_scores[:, 0], data_mod_scores[:, 1], color = 'red', label = 'data')
ax.plot(shuffled_mod_scores[:, 0], shuffled_mod_scores[:, 1], color = 'blue', label = 'shuffled')
ax.set_xlabel('gamma')
ax.set_ylabel('modularity score')
ax.set_title(name)
ax.legend()
fig.savefig(output_dir+name+".png", format="png")

[[community], mod_score] = bct.modularity_louvain_dir(r, gamma=data_mod_scores[max_gamma][0], hierarchy=True)

num_communities = len(np.unique(community))
list_communities = []
for i in range(num_communities):
    list_communities.append([])

for i in range(len(cortical_areas)):
    community_number = community[i]
    list_communities[community_number-1].append(cortical_areas[i])
for module in list_communities:
    print(module)