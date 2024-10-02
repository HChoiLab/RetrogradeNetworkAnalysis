#!/usr/bin/env python
# coding: utf-8


import logging
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

"""
Published on Tues Oct 1 11:22:46 2024

@author: Eric Li, Hannah Choi
"""
"""
This code takes previously generated hierarchy scores & global hierarchy score data from the "Results" folder, and makes figures.
"""


# # Import results

input_dir = r'./results/' 
output_dir = r'./Output/'
xls1=pd.ExcelFile(input_dir+"hierarchy_summary.xlsx")
xls2=pd.ExcelFile(input_dir+"gh_comparison_shuffled_CreConf.xlsx")

# # Hierarchy of all cortical + thalamic + basal ganglia regions

df=pd.read_excel(xls1)
df = df.sort_values('CC+TCCT+CB iterated')
areas = df['areas']

areas_isocortex = df[df.CortexThalamusBG == 'C']['areas']
areas_T = df[df.CortexThalamusBG != 'B']['areas']

areas_CC = df['areas'].copy()
mask = np.in1d(areas_CC, areas_isocortex)
areas_CC.loc[~mask] = np.nan
areas_vec = np.arange(0,len(areas))
areas_vec_CC = np.arange(0,len(areas))
areas_vec_CC = [x for i, x in enumerate(areas_vec_CC) if mask[i]]

areas_TCCT = df['areas'].copy()
mask = np.in1d(areas_TCCT, areas_T)
areas_TCCT.loc[~mask] = np.nan
areas_vec = np.arange(0,len(areas))
areas_vec_TCCT = np.arange(0,len(areas))
areas_vec_TCCT = [x for i, x in enumerate(areas_vec_TCCT) if mask[i]]

hs = df['CC+TCCT+CB iterated']
hs_CC = df[df.CortexThalamusBG == 'C']['CC iterated']
hs_TCCT = df[df.CortexThalamusBG != 'B']['CC+TCCT iterated']

f = plt.figure()
plt.plot(hs_CC, areas_vec_CC, 'ro', label = 'CC')
plt.plot(hs_TCCT, areas_vec_TCCT, 'go', label = 'CC+TCCT')
plt.plot(hs,areas_vec,'bo', label = 'CC+TCCT+CB')
plt.xlabel('hierarchy score')
plt.ylabel('areas')
plt.xlim(-2.5,1.5)
plt.yticks(areas_vec, areas)
plt.title('CC+TCCT+CB hierarchy')
f.set_size_inches(6, 10.5)
plt.legend(loc='lower right')
# f.savefig(output_dir+"hierarchy.svg", format="svg")
plt.show() 


# # Intra-module hierarchy of visual module

df=pd.read_excel(xls1,'hierarchy_visual')
df = df.sort_values('CC+TCCT iterated')

areas = df['areas']
areas_vec = np.arange(0,len(areas))
hs = df['CC+TCCT iterated']

f = plt.figure()
plt.plot(hs,areas_vec,'ro')
plt.xlabel('hierarchy score')
plt.ylabel('areas')
plt.xlim(-1.1,1)
plt.yticks(areas_vec, areas)
plt.title('Visual module hierarchy')
f.set_size_inches(4, 6)
#f.savefig(output_dir+"visual.png", format="png")
plt.show()


# # Inter-module cortical hierarchy

''' Inter-module cortical hierarchy'''

df=pd.read_excel(xls1,'hierarchy_intermodule')
df = df.sort_values('CC+TCCT iterated')

areas = df['areas']
areas_vec = np.arange(0,len(areas))
hs = df['CC+TCCT iterated']

f = plt.figure()
plt.plot(hs,areas_vec,'ro')
plt.xlabel('hierarchy score')
plt.ylabel('areas')
plt.xlim(-1.1,1)
plt.yticks(areas_vec, areas)
plt.title('Inter-module hierarchy')
f.set_size_inches(4, 6)
#f.savefig(output_dir+"inter-module.png", format="png")
plt.show()


# # Global hierarchy scores: comparison to shuffled data

df_CC = pd.read_excel(xls2,'CC')
CC_shuffled = df_CC['iter shuffled hg (cortex)']
CC_data = df_CC['iter data hg (cortex)'][0]
zscore_CC = (CC_data-np.mean(CC_shuffled))/np.std(CC_shuffled) 

df_TCCT = pd.read_excel(xls2,'CC+TCCT')
TCCT_shuffled = df_TCCT['iter shuffled hg (cortex+thal)']
TCCT_data = df_TCCT['iter data hg (cortex+thal)'][0]
zscore_TCCT = (TCCT_data-np.mean(TCCT_shuffled))/np.std(TCCT_shuffled) 

df_CB = pd.read_excel(xls2,'CC+TCCT+CB')
CB_shuffled = df_CB['iter shuffled hg (cortex+thal+bg)']
CB_data = df_CB['iter data hg (cortex+thal+bg)'][0]
zscore_CB = (CB_data-np.mean(CB_shuffled))/np.std(CB_shuffled) 

fig,ax=plt.subplots()
ax.hist(CC_shuffled, bins=10, label='CC',color='red')
ax.axvline(x=CC_data,linestyle='--',color='red')
ax.hist(TCCT_shuffled, bins=25, label='CC+CTTC')
ax.axvline(x=TCCT_data,linestyle='--')
ax.hist(CB_shuffled, bins=25, label='CC+CTTC+CB', color = 'orange')
ax.axvline(x=CB_data,linestyle='--', color='orange')
ax.set_xlabel('global hierarchy score (antero + retro)',fontsize=16)
ax.set_ylabel('counts',fontsize=16)
ax.set_title('z-score(CC)={:.2f}, z-score(CC+CTTC)={:.2f}, z-score(CC+CTTC+CB)={:.2f}'.format(zscore_CC, zscore_TCCT, zscore_CB), fontsize='small')
ax.legend(loc='upper right')
# fig.savefig(output_dir+"shuffled.svg", format="svg")