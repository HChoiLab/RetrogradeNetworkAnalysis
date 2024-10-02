#!/usr/bin/env python
# coding: utf-8

import logging
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from IterativeMethod import iterativeCB

"""
Published on Tues Oct 1 16:40:37 2024

@author: Eric Li, Hannah Choi
"""
"""
This code finds hierarchy scores of cortical, thalamic, and basal ganglia areas via iterations, based on
cortico-basal ganglia (CB) connections, in addition to cortico-thalamic network.
"""

# Set input and output directories

input_dir = r'../Input/'     # Directory with the file "CB_clusters.csv" 
output_dir = r'../Output/'   # Directory to save the outputs from the experimental data

CreConf = 1                 # 1 if using CC hierarchy with Cre-confidence; 0 if not

# Read the excel file with source-target-creline pairs and their cluster numbers.

df=pd.read_csv(input_dir+"CB_clusters.csv")

dfV = df[['source_area','target','Cre','tree']]
dfV = dfV.reset_index(drop=True)
dfV.columns = ['source', 'target', 'creline', 'clu']
dfV['clu'] = dfV['clu'].apply(np.int64)

num_clu = len(dfV["clu"].unique())

dfVB = dfV.reset_index(drop=True)

# Set mappings to all feedforward.

jmax = 2**num_clu-1

# Define functions needed

c0 = 2**num_clu

def ffb (cls):
    """Direction of each CB cluster"""
    b=(bin(c0+jmax)[-(cls)])
    return (2*int(b)-1)
    
def hrf (area):
    '''Hierarchy score of each area'''
    return np.mean(dfVB[dfVB.target == area].ffb)

# Produce expanded data frame with  FF/FB, hierarchy values as source & target for each pair of CB connections

dfVB["ffb"] = dfVB["clu"].apply(ffb)
dfVB["hr_t"]=dfVB["target"].apply(hrf)

dfVB.to_excel(output_dir+'inputexpanded_CB.xlsx')

# Find hierarchy score for each of basal ganglia areas

areas = dfVB["target"].unique()
n_areas=len(areas) 
hr=list(range(0,n_areas))

for i in range(0,n_areas):
    hr[i]=np.mean(dfVB[dfVB.target == areas[i]].ffb)
    
data=[areas,hr]
data=np.transpose(data)
columns = ['areas','h']
dfiB = pd.DataFrame(data,columns=columns)
dfiB.head()

# Iterate basal ganglia + thalamic + cortical hierarchy scores

n_iter = 20

dfiCT = pd.read_excel(output_dir+"TCCT_CCconf_iter.xlsx")
dfiCT['h'] = dfiCT[n_iter]

dfVC = pd.read_excel(output_dir+"inputexpanded_CC.xlsx")
dfVT = pd.read_excel(output_dir+"inputexpanded_TCCT.xlsx")
dfVCT = pd.concat([dfVC, dfVT]).reset_index(drop=True)
dfVCT = dfVCT[["source","target","ffb_nc","ffb_c"]]

if CreConf == 0:
    dfVCT["ffb"] = dfVCT["ffb_nc"]
elif CreConf == 1:
    dfVCT["ffb"] = dfVCT["ffb_c"]

dfiB = dfiB[["areas","h"]]
dfiCT = dfiCT[["areas","h"]]
dfVB = dfVB[["source","target","ffb"]]
dfVCT = dfVCT[["source","target","ffb"]]

'''Iterations ''' 
hr_iter = iterativeCB(dfiCT, dfVCT, dfiB, dfVB, n_iter)
iteration=np.arange(0,n_iter+1,1)
n_area=np.shape(hr_iter)[0]
allareas = hr_iter["areas"].unique()

'''Figure of hierarchy score iterations'''

fig,ax=plt.subplots()
for i in range(0,n_area):
    y=np.squeeze(np.asarray(hr_iter[hr_iter.areas==allareas[i]].iloc[:,1::1]))
    ax.plot(iteration,y)
ax.set_xlim([0, n_iter])
ax.set_xticks(np.arange(0, n_iter, step=5))
ax.set_xlabel('iter')
ax.set_ylabel('hierarchy value')
plt.show()

'''Figure showing correlation between hierarchy scores before & after iterations'''

hr_final = hr_iter[:][n_iter]
hr_initial = hr_iter[:][0]   
f = plt.figure()
plt.plot(hr_initial,hr_final,'ro')
plt.xlabel('initial hierarchy (conf)')
plt.ylabel('final hierarchy (conf)')
plt.title('r='+str(np.corrcoef(hr_initial, hr_final)[0, 1]))
plt.show()

'''Save hierarchy scores before and after iteration'''
if CreConf == 0:
    dfC = pd.read_excel(output_dir+"CC_noconf_iter.xlsx")
elif CreConf == 1:
    dfC = pd.read_excel(output_dir+"CC_conf_iter.xlsx")

for i_area in range(0,n_area):
    if hr_iter['areas'][i_area] in list(dfC['areas']):
        hr_iter.loc[i_area,'CortexThalamusBG'] = 'C'
    elif hr_iter['areas'][i_area] in list(dfiCT['areas']):
        hr_iter.loc[i_area,'CortexThalamusBG'] = 'T'
    else:
        hr_iter.loc[i_area,'CortexThalamusBG'] = 'B'

hr_iter = hr_iter[['areas','CortexThalamusBG', 0,n_iter]]  

if CreConf == 0:
    hr_iter.to_excel(output_dir+'TCCTCBnoconf_iter.xlsx')
elif CreConf == 1:
    hr_iter.to_excel(output_dir+'TCCTCBconf_iter.xlsx')

# Print out global hierarchy scores for the CC+TCCT+CB connectivity data before and after iteration

dfi_CTB = hr_iter[["CortexThalamusBG","areas",0,n_iter]]
dfV_CB = dfVB[["source","target","ffb"]]
dfV_CT = dfVCT[['source','target','ffb']]
dfi_CTB = dfi_CTB.rename(columns={0: "h0", n_iter:"h_iter"})

dfi_cortex_thalamus1 = dfi_CTB[(dfi_CTB.CortexThalamusBG == 'C') | (dfi_CTB.CortexThalamusBG == 'T')]
dfi_cortex_thalamus1 = dfi_cortex_thalamus1[['areas','h_iter']]
dfV_CT = dfV_CT.join(dfi_cortex_thalamus1.set_index('areas'), on ='source')
dfV_CT=dfV_CT.rename(columns={"h_iter": "hs"})
dfV_CT = dfV_CT.join(dfi_cortex_thalamus1.set_index('areas'), on ='target')
dfV_CT=dfV_CT.rename(columns={"h_iter": "ht"})
dfV_CT = dfV_CT.dropna() 
hg_CT_1 = dfV_CT.ffb*(dfV_CT.ht - dfV_CT.hs)

dfi_bg1=dfi_CTB[(dfi_CTB.CortexThalamusBG == 'B')]
dfi_bg1 = dfi_bg1[['areas','h_iter']]
dfV_CB = dfV_CB.join(dfi_bg1.set_index('areas'), on ='target')
dfV_CB=dfV_CB.rename(columns={"h_iter": "ht"})
dfV_CB = dfV_CB.join(dfi_cortex_thalamus1.set_index('areas'), on ='source')
dfV_CB=dfV_CB.rename(columns={"h_iter": "hs"})
dfV_CB = dfV_CB.dropna() 
hg_CB_1 = dfV_CB.ffb*(dfV_CB.ht- dfV_CB.hs)

hg_cortexthalamus_CB_iter = np.mean(hg_CT_1)
hg_CB_iter = np.mean(hg_CT_1.append(hg_CB_1))

dfV_CB = dfVB[["source","target","ffb"]]
dfV_CT = dfVCT[['source','target','ffb']]

dfi_cortex_thalamus1 = dfi_CTB[(dfi_CTB.CortexThalamusBG == 'C') | (dfi_CTB.CortexThalamusBG == 'T')]
dfi_cortex_thalamus1 = dfi_cortex_thalamus1[['areas','h0']]
dfV_CT = dfV_CT.join(dfi_cortex_thalamus1.set_index('areas'), on ='source')
dfV_CT=dfV_CT.rename(columns={"h0": "hs"})
dfV_CT = dfV_CT.join(dfi_cortex_thalamus1.set_index('areas'), on ='target')
dfV_CT=dfV_CT.rename(columns={"h0": "ht"})
dfV_CT = dfV_CT.dropna() 
hg_CT_1 = dfV_CT.ffb*(dfV_CT.ht- dfV_CT.hs)

dfi_bg1=dfi_CTB[(dfi_CTB.CortexThalamusBG == 'B')]
dfi_bg1 = dfi_bg1[['areas','h0']]
dfV_CB = dfV_CB.join(dfi_bg1.set_index('areas'), on ='target')
dfV_CB=dfV_CB.rename(columns={"h0": "ht"})
dfV_CB = dfV_CB.join(dfi_cortex_thalamus1.set_index('areas'), on ='source')
dfV_CB=dfV_CB.rename(columns={"h0": "hs"})
dfV_CB = dfV_CB.dropna() 
hg_CB_1 = dfV_CB.ffb*(dfV_CB.ht- dfV_CB.hs)

hg_bg_CB_init = np.mean(hg_CB_1)
hg_cortexthalamus_CB = np.mean(hg_CT_1)
hg_CB_init = np.mean(hg_CT_1.append(hg_CB_1))

print('hg of CB BG='+str(hg_bg_CB_init))
print('hg of CC+TCCT+CB before iterations ='+str(hg_CB_init)) 
print('hg of CC+TCCT+CB iterate cortex+thalamus ='+str(hg_cortexthalamus_CB_iter))
print('hg of CC+TCCT+CB iterate cortex+thalamus+bg='+str(hg_CB_iter))


'''Save global hierarchy scores'''
newDF= pd.DataFrame([])
newDF=newDF.append(pd.DataFrame({'hg_bg_CB_init':hg_bg_CB_init, 'hg_CB_init':hg_CB_init,
                                 'hg_cortexthalamus_CB_iter':hg_cortexthalamus_CB_iter, 
                                 'hg_CB_iter':hg_CB_iter},index=[0]))
newDF.to_excel(output_dir+'ghs_CB.xlsx')