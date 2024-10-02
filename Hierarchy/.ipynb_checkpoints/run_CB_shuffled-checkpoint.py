#!/usr/bin/env python
# coding: utf-8

import logging
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from IterativeMethod import iterativeCB
"""
Published on Tues October 1 11:03:14 2024

@author: Eric Li, Hannah Choi
"""
"""
This code generates shuffled versions of the CB connectivity data and 
computes the global hierarchy scores of the shuffled data.
"""

"""
Note: Run run_TCCT_shuffled.py beforehand!
"""

# Set input and output directories 

input_dir = r'../Input/'        # Directory with the file "CB_clusters.csv"
input_dir2 = r'../Output/'               # Directory with the file "ghc_CB.xls"
output_dir = r'../Output/shuffled/'   # Directory to save the ouputs from the shuffled experimental data

CreConf = 1                 # 1 if using CC hierarchy with Cre-confidence; 0 if not

# Read the excel file with CB source-target-creline pairs and their cluster numbers.

df=pd.read_csv(input_dir+"CB_clusters.csv")

dfV = df[['source_area','target','Cre','tree']]
dfV = dfV.reset_index(drop=True)
dfV.columns = ['source', 'target', 'creline', 'clu']
dfV['clu'] = dfV['clu'].apply(np.int64)

num_clu = len(dfV["clu"].unique())

dfVCB = dfV.reset_index(drop=True)

# Find global hierarchy scores of shuffled CB connectivity data

n_iter = 10
n_shuffle = 100
by_creline = 0
line_list= dfVCB["creline"].unique()

hr_cb_shuffled = np.zeros(n_shuffle)
hr_iter_ct_shuffled = np.zeros(n_shuffle)
hr_iter_cb_shuffled = np.zeros(n_shuffle)

conf_shuffled = np.zeros(n_shuffle)

for i_shuffle in range(0,n_shuffle):
    
    print('i_shuffle='+str(i_shuffle))
    
    dfVB_shuffled = dfVCB
    
    if by_creline == 0:        
        source_list = dfVB_shuffled.source
        target_list = dfVB_shuffled.target
        source_shuffled = source_list.sample(frac=1).reset_index(drop=True)
        target_shuffled = target_list.sample(frac=1).reset_index(drop=True)
        source_shuffled.index = source_list.index
        target_shuffled.index = target_list.index 
        dfVB_shuffled.loc[source_list.index, "source"]=np.array(source_shuffled)
        dfVB_shuffled.loc[target_list.index, "target"]=np.array(target_shuffled)
        
    elif by_creline == 1:
        for i in range(0,len(line_list)):
            source_list = dfVB_shuffled[(dfVB_shuffled.creline == str(line_list[i]))].source
            target_list = dfVB_shuffled[(dfVB_shuffled.creline == str(line_list[i]))].target
            source_shuffled = source_list.sample(frac=1).reset_index(drop=True)
            target_shuffled = target_list.sample(frac=1).reset_index(drop=True)
            source_shuffled.index = source_list.index
            target_shuffled.index = target_list.index        
            dfVB_shuffled.loc[source_list.index, "source"]=np.array(source_shuffled)
            dfVB_shuffled.loc[target_list.index, "target"]=np.array(target_shuffled)

    dfVB_shuffled['clu'] = dfVB_shuffled['clu'].astype('Int64')

    source_areas = dfVB_shuffled["source"].unique()
    target_areas = dfVB_shuffled["target"].unique()
    
    areas = source_areas
    n_areas=len(areas)

    n = num_clu
        
    c0 = 2**n
    
    """Define functions needed"""
    
    def ff_or_fb (cls):
        """Direction of each CB cluster"""
        return 1
        
    def hrf (area):
        '''Hierarchy score of each area'''
        return np.mean(dfVB_shuffled[dfVB_shuffled.target == area].ffb)

###########################################################################
    
    """Produce an expanded data frame with FF/FB, hierarchy values as source & target 
    for each pair of CB connections"""
    
    dfVB_shuffled["ffb"]=dfVB_shuffled["clu"].apply(ff_or_fb)
    dfVB_shuffled.loc[:,"hr_t"]=dfVB_shuffled["target"].apply(hrf)

    ###########################################################################
    '''Finding initial hierarchy score of each BG area'''    

    '''Iterate bg + thalamic + cortical hierarchy scores'''
    areas = target_areas
    n_areas=len(areas) 
    hr=list(range(0,n_areas))
    
    for i in range(0,n_areas):
        hr[i]=np.mean(dfVB_shuffled[dfVB_shuffled.target == areas[i]].ffb)

    data=[areas,hr]
    data=np.transpose(data)
    columns = ['areas','h']
    dfiB = pd.DataFrame(data,columns=columns)
    ##########################################################
    '''Load results from CC+TCCT hierarchy'''

    dfiCT = pd.read_excel(output_dir+'TCCT_CCconf_iter'+str(i_shuffle)+'.xlsx')  
    dfiCT['h'] = dfiCT[n_iter]
    
    dfVC = pd.read_excel(output_dir+'inputexpanded_CC_shuffled'+str(i_shuffle)+'.xlsx')
    dfVT = pd.read_excel(output_dir+'inputexpanded_TCCT_shuffled'+str(i_shuffle)+'.xlsx')
    dfVCT = pd.concat([dfVC, dfVT]).reset_index(drop=True)
    dfVCT = dfVCT[["source","target","ffb_nc","ffb_c"]]

    if CreConf == 0:
        dfVCT["ffb"] = dfVCT["ffb_nc"]
    elif CreConf == 1:
        dfVCT["ffb"] = dfVCT["ffb_c"]
    ##########################################################    
    
    dfiB = dfiB[["areas","h"]]
    dfiCT = dfiCT[["areas","h"]]
    dfVB = dfVB_shuffled[["source","target","ffb"]]
    dfVCT = dfVCT[["source","target","ffb"]]
    
    hr_iter = iterativeCB(dfiCT, dfVCT, dfiB, dfVB, n_iter)    
    
    iteration=np.arange(0,n_iter+1,1)
    n_area=np.shape(hr_iter)[0]
    allareas = hr_iter["areas"].unique()

    if CreConf == 0:
        dfC = pd.read_excel(output_dir+"CCshuffled_noconf_iter"+str(i_shuffle)+".xlsx")
    elif CreConf == 1:
        dfC = pd.read_excel(output_dir+"CCshuffled_conf_iter"+str(i_shuffle)+".xlsx")
    
    for i_area in range(0,n_area):
        if hr_iter['areas'][i_area] in list(dfC['areas']):
            hr_iter.loc[i_area,'CortexThalamusBG'] = 'C'
        elif hr_iter['areas'][i_area] in list(dfiCT['areas']):
            hr_iter.loc[i_area,'CortexThalamusBG'] = 'T'
        else:
            hr_iter.loc[i_area,'CortexThalamusBG'] = 'B'
    
    hr_iter = hr_iter[['areas','CortexThalamusBG', 0,n_iter] ]  
   
    if CreConf == 1:
        hr_iter.to_excel(output_dir+'CB_TCCT_CCconf_iter'+str(i_shuffle)+'.xlsx')
    elif CreConf == 0:
        hr_iter.to_excel(output_dir+'CB_TCCT_CCnoconf_iter'+str(i_shuffle)+'.xlsx')
    ###########################################################################


    ###########################################################################
    '''global hierarchy score of the shuffled data before & after iteration'''
    
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
    
    hr_cb_shuffled[i_shuffle] = hg_CB_init
    hr_iter_ct_shuffled[i_shuffle] = hg_cortexthalamus_CB_iter
    hr_iter_cb_shuffled[i_shuffle] = hg_CB_iter

pd.DataFrame(hr_cb_shuffled).to_excel(output_dir+'shuffled_hg_CB_init.xlsx')
pd.DataFrame(hr_iter_ct_shuffled).to_excel(output_dir+'shuffled_hg_CB_iter_cortexthalamus.xlsx')
pd.DataFrame(hr_iter_cb_shuffled).to_excel(output_dir+'shuffled_hg_CB_iter.xlsx')

# Plot global hierarchy scores of 100 shuffled data with the global hierarchy score of the original data

"""Global hierarchy scores of the original CB connectivity"""

df_hg_CB = pd.read_excel(input_dir2+'ghs_CB.xlsx')

hg_all_init = df_hg_CB["hg_CB_init"][0]
hg_cortex_iter = df_hg_CB["hg_cortexthalamus_CB_iter"][0]
hg_all_iter = df_hg_CB["hg_CB_iter"][0] 

hm1 = (hg_all_init-np.mean(hr_cb_shuffled))/np.std(hr_cb_shuffled)              # Z-score for thalamus+cortex+bg before iteration
hm2 = (hg_cortex_iter-np.mean(hr_iter_cb_shuffled))/np.std(hr_iter_cb_shuffled)   # Z-score for thalamus+cortex after iteration
hm3 = (hg_all_iter-np.mean(hr_iter_cb_shuffled))/np.std(hr_iter_cb_shuffled)    # Z-score for thalamus+cortex+bg after iteration

""" Figure showing global hierarchy scores of shuffled data & original CB data before iteration, for cortex+thalamus+bg """
fig,ax=plt.subplots()
bins=25
ax.hist(hr_cb_shuffled, bins=bins)
ax.axvline(x=hg_all_init,linestyle='--')
ax.set_xlabel('global hierarchy score',fontsize=16)
ax.set_ylabel('counts',fontsize=16)
ax.set_title('hm='+str(hm1))

""" Figure showing global hierarchy scores of shuffled data & original CB data after iteration, for cortex+thalamus only """
fig,ax=plt.subplots()
bins=25
ax.hist(hr_iter_cb_shuffled, bins=bins)
ax.axvline(x=hg_cortex_iter,linestyle='--')
ax.set_xlabel('global hierarchy score',fontsize=16)
ax.set_ylabel('counts',fontsize=16)
ax.set_title('hm='+str(hm2))

""" Figure showing global hierarchy scores of shuffled data & original CB data after iteration, for cortex+thalamus+bg """
fig,ax=plt.subplots()
bins=25
ax.hist(hr_iter_cb_shuffled, bins=bins)
ax.axvline(x=hg_all_iter,linestyle='--')
ax.set_xlabel('global hierarchy score',fontsize=16)
ax.set_ylabel('counts',fontsize=16)
ax.set_title('hm='+str(hm3))