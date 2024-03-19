#!/usr/bin/env python
# coding: utf-8


import logging
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from func_unsupervised_CC_retro import fit_retro
from func_unsupervised_CC_antero import fit_antero
from IterativeMethod import iterativeCC

"""
Published on Sat Mar 2 17:23:54 2024

@author: Eric Li, Hannah Choi
"""
"""
This code generates shuffled versions of the CC connectivity data and 
computes the global hierarchy scores of the shuffled data.
"""


# # Set input and output directories


input_dir = r'./Input/'                     # Directory with the file "AnteroRetro_CC_TC_CT_clusters.xlsx"
input_dir2 = r'./Output/'                   # Directory with the file "ghs_CC.xlsx"
output_dir = r'./Output/shuffled/'          # Directory to save the ouputs from the shuffled experimental data


# # Read the excel file with source-target-creline pairs and their cluster numbers.


xls=pd.ExcelFile(input_dir+"AnteroRetro_CC_TC_CT_clusters.xlsx")
df=pd.read_excel(xls,'anteroCCTC_retroCCCT')

df_antero=df[(df.Antero_Retro == "A")] 
df_antero=df_antero[(df_antero.hemi == "ipsi")&(df_antero.creline != "C57BL/6J / Emx1")&(df_antero.target != "VISC")&(df_antero.source != "SSp-un")&(df_antero.target != "SSp-un")&(df_antero.source != "VISC")] # Consider ipsilateral, Cre connectivity data only
df_antero = df_antero[(df_antero["Target Major Division"] == "isocortex")&(df_antero["Source Major Division"] == "isocortex")] # Consider C-C connections only

dfV1 = df_antero[['Antero_Retro','source','target','creline','AnteroCluster']]
dfV1 = dfV1.reset_index(drop=True)
dfV1['AnteroCluster'] = dfV1['AnteroCluster'].apply(np.int64)

source_areas1 = dfV1["source"].unique()
target_areas1 = dfV1["target"].unique()

df_retro=df[(df.Antero_Retro == "R")] 
df_retro = df_retro[(df_retro["Target Major Division"] == "isocortex")&(df_retro["Source Major Division"] == "isocortex")]  # Consider C-C connections only

dfV2 = df_retro[['Antero_Retro','source','target','creline','RetroCluster']]
dfV2 = dfV2.reset_index(drop=True)
dfV2['RetroCluster'] = dfV2['RetroCluster'].apply(np.int64)

source_areas2 = dfV2["source"].unique()
target_areas2 = dfV2["target"].unique()  

dfV_antero = dfV1[["source","target","creline","AnteroCluster", "Antero_Retro"]].copy()
dfV_retro = dfV2[["source","target","creline","RetroCluster", "Antero_Retro"]].copy()


# # Find global hierarchy scores of shuffled CC connectivity data


n_iter = 10 # or 5, for fast process
n_shuffle = 100
by_creline = 1   # 1 if shuffle within each Cre-line; 0 if shuffle across Cre-lines 
line_list_antero = dfV_antero["creline"].unique()
line_list_retro = dfV_retro["creline"].unique()

hr_init_shuffled = np.zeros(n_shuffle)
hr_iter_shuffled = np.zeros(n_shuffle)
hrc_iter_shuffled = np.zeros(n_shuffle)
hrc_init_shuffled = np.zeros(n_shuffle)

for i_shuffle in range(0,n_shuffle):
    dfV_shuffled_antero = dfV_antero
    dfV_shuffled_retro = dfV_retro
    
    if by_creline == 0:
        source_list_antero = dfV_shuffled_antero.source
        target_list_antero = dfV_shuffled_antero.target
        source_shuffled_antero = source_list_antero.sample(frac=1).reset_index(drop=True)
        target_shuffled_antero = target_list_antero.sample(frac=1).reset_index(drop=True)
        source_shuffled_antero.index = source_list_antero.index
        target_shuffled_antero.index = target_list_antero.index 
        dfV_shuffled_antero.loc[source_list_antero.index, "source"]=np.array(source_shuffled_antero)
        dfV_shuffled_antero.loc[target_list_antero.index, "target"]=np.array(target_shuffled_antero)
        
        source_list_retro = dfV_shuffled_retro.source
        target_list_retro = dfV_shuffled_retro.target
        source_shuffled_retro = source_list_retro.sample(frac=1).reset_index(drop=True)
        target_shuffled_retro = target_list_retro.sample(frac=1).reset_index(drop=True)
        source_shuffled_retro.index = source_list_retro.index
        target_shuffled_retro.index = target_list_retro.index 
        dfV_shuffled_retro.loc[source_list_retro.index, "source"]=np.array(source_shuffled_retro)
        dfV_shuffled_retro.loc[target_list_retro.index, "target"]=np.array(target_shuffled_retro)
        
    elif by_creline == 1:
        for i in range(0,len(line_list_antero)):
            source_list_antero = dfV_shuffled_antero[(dfV_shuffled_antero.creline == str(line_list_antero[i]))].source
            target_list_antero = dfV_shuffled_antero[(dfV_shuffled_antero.creline == str(line_list_antero[i]))].target
            source_shuffled_antero = source_list_antero.sample(frac=1).reset_index(drop=True)
            target_shuffled_antero = target_list_antero.sample(frac=1).reset_index(drop=True)
            source_shuffled_antero.index = source_list_antero.index
            target_shuffled_antero.index = target_list_antero.index        
            dfV_shuffled_antero.loc[source_list_antero.index, "source"]=np.array(source_shuffled_antero)
            dfV_shuffled_antero.loc[target_list_antero.index, "target"]=np.array(target_shuffled_antero)
            
        for i in range(0,len(line_list_retro)):
            source_list_retro = dfV_shuffled_retro[(dfV_shuffled_retro.creline == str(line_list_retro[i]))].source
            target_list_retro = dfV_shuffled_retro[(dfV_shuffled_retro.creline == str(line_list_retro[i]))].target
            source_shuffled_retro = source_list_retro.sample(frac=1).reset_index(drop=True)
            target_shuffled_retro = target_list_retro.sample(frac=1).reset_index(drop=True)
            source_shuffled_retro.index = source_list_retro.index
            target_shuffled_retro.index = target_list_retro.index        
            dfV_shuffled_retro.loc[source_list_retro.index, "source"]=np.array(source_shuffled_retro)
            dfV_shuffled_retro.loc[target_list_retro.index, "target"]=np.array(target_shuffled_retro)

    frames = [dfV_shuffled_antero, dfV_shuffled_retro]
    dfV_shuffled = pd.concat(frames)
    dfV_shuffled['AnteroCluster'] = dfV_shuffled['AnteroCluster'].astype('Int64')
    dfV_shuffled['RetroCluster'] = dfV_shuffled['RetroCluster'].astype('Int64')

    source_areas = dfV_shuffled["source"].unique()
    target_areas = dfV_shuffled["target"].unique()
    
    source_areas1 = set(source_areas)
    target_areas1 = set(target_areas)

    areas = np.array(list(source_areas1.intersection(target_areas1)))
    n_areas=len(areas)
    
    num_antero_clu = len(dfV_shuffled["AnteroCluster"].unique())-1
    num_retro_clu = len(dfV_shuffled["RetroCluster"].unique())-1
    dfV_shuffled = dfV_shuffled.reset_index(drop=True)
        
    hierarchy_vals_antero = fit_antero(dfV_shuffled_antero)
        
    jmax_raw_antero, jmax_antero = np.argmax(hierarchy_vals_antero, axis=0)
    jmax_raw_val_antero = hierarchy_vals_antero[jmax_raw_antero][0]
    jmax_val_antero = hierarchy_vals_antero[jmax_antero][1]
    logging.debug("RESULTS")
    n_a = num_antero_clu
    logging.debug("(jmax_raw_antero, val) = ({:0{n}b}, {:.3f})".format(jmax_raw_antero, jmax_raw_val_antero, n=n_a))
    logging.debug("(jmax_antero,     val) = ({:0{n}b}, {:.3f})".format(jmax_antero, jmax_val_antero, n=n_a))

    results_antero = dict(jmax_antero=bin(2**n_a+jmax_antero),
               jmax_val_antero=jmax_val_antero,
               jmax_raw_antero=bin(2**n_a+jmax_raw_antero),
               jmax_raw_val_antero=jmax_raw_val_antero)
    
    hierarchy_vals_retro = fit_retro(dfV_shuffled_retro)
        
    jmax_raw_retro, jmax_retro = np.argmax(hierarchy_vals_retro, axis=0)
    jmax_raw_val_retro = hierarchy_vals_retro[jmax_raw_retro][0]
    jmax_val_retro = hierarchy_vals_retro[jmax_retro][1]
    logging.debug("RESULTS")
    n_r = num_retro_clu
    logging.debug("(jmax_raw_retro, val) = ({:0{n}b}, {:.3f})".format(jmax_raw_retro, jmax_raw_val_retro, n=n_r))
    logging.debug("(jmax_retro,     val) = ({:0{n}b}, {:.3f})".format(jmax_retro, jmax_val_retro, n=n_r))

    results_retro = dict(jmax_retro=bin(2**n_r+jmax_retro),
               jmax_val_retro=jmax_val_retro,
               jmax_raw_retro=bin(2**n_r+jmax_raw_retro),
               jmax_raw_val_retro=jmax_raw_val_retro)
    
    print('i_shuffle='+str(i_shuffle))
    
    ###########################################################################    
    """Define functions needed"""
    
    c0r=2**(num_retro_clu)
    c0a=2**(num_antero_clu)
    
    def ff_or_fb_retro (cls):
        """Direction of each retrograde CC cluster with Cre-confidence"""
        b=(bin(c0r+jmax_retro)[-(cls)])
        return (2*int(b)-1)
    
    def ff_or_fb_retro_nc (cls):
        """Direction of each retrograde CC cluster without Cre-confidence"""
        b=(bin(c0r+jmax_raw_retro)[-(cls)])
        return (2*int(b)-1)
    
    def ff_or_fb_antero(cls):
        """Direction of each anterograde CC cluster with Cre-confidence"""
        b=(bin(c0a+jmax_antero)[-(cls)])
        return (2*int(b)-1)

    def ff_or_fb_antero_nc(cls):
        """Direction of each anterograde CC cluster without Cre-confidence"""
        b=(bin(c0a+jmax_raw_antero)[-(cls)])
        return (2*int(b)-1)
    
    def cre_confidence1(df):
        """Returns confidence of cre lines"""
        func = lambda x: 1 - np.abs(x.mean())
        return df.groupby('creline')['ffb_c'].transform(func)
    
    def hrf (area):
        '''Hierarchy score of each area without Cre-confidence'''
        return ((-np.mean(dfV_shuffled[dfV_shuffled.source == area].ffb_nc)
                 +np.mean(dfV_shuffled[dfV_shuffled.target == area].ffb_nc))/2)
    
    def hrcf (area):
        '''Hierarchy score of each area with Cre-confidence'''
        return ((-np.mean(dfV_shuffled[dfV_shuffled.source == area].ffb_c*dfV_shuffled[dfV_shuffled.source == area].conf)
                 +np.mean(dfV_shuffled[dfV_shuffled.target == area].ffb_c*dfV_shuffled[dfV_shuffled.target == area].conf))/2)
    ###########################################################################
    
    ###########################################################################
    """Produce expanded data frame with  FF/FB, Cre-confidence, hierarchy values as source & target 
    for each pair of CC connections"""
    
    dfV_shuffled.loc[(dfV_shuffled.Antero_Retro == "R"),"ffb_c"]=dfV_shuffled[(dfV_shuffled.Antero_Retro == "R")].RetroCluster.apply(ff_or_fb_retro)
    dfV_shuffled.loc[(dfV_shuffled.Antero_Retro == "R"),"ffb_nc"]=dfV_shuffled[(dfV_shuffled.Antero_Retro == "R")].RetroCluster.apply(ff_or_fb_retro_nc)
    dfV_shuffled.loc[(dfV_shuffled.Antero_Retro == "A"),"ffb_c"]=dfV_shuffled[(dfV_shuffled.Antero_Retro == "A")].AnteroCluster.apply(ff_or_fb_antero)
    dfV_shuffled.loc[(dfV_shuffled.Antero_Retro == "A"),"ffb_nc"]=dfV_shuffled[(dfV_shuffled.Antero_Retro == "A")].AnteroCluster.apply(ff_or_fb_antero_nc)
    dfV_shuffled.loc[:, "conf"] = cre_confidence1(dfV_shuffled)
    dfV_shuffled.loc[:,"hrc_s"]=dfV_shuffled["source"].apply(hrcf)
    dfV_shuffled.loc[:,"hrc_t"]=dfV_shuffled["target"].apply(hrcf)
    dfV_shuffled.loc[:,"hr_s"]=dfV_shuffled["source"].apply(hrf)
    dfV_shuffled.loc[:,"hr_t"]=dfV_shuffled["target"].apply(hrf)
    
    dfV_shuffled.to_excel(output_dir+'inputexpanded_CC_shuffled'+str(i_shuffle)+'.xlsx')
    ###########################################################################
    
    hrs=list(range(0,n_areas))
    hrt=list(range(0,n_areas))
    hr=list(range(0,n_areas))
    hrc=list(range(0,n_areas))
    for i in range(0,n_areas):
        hrs[i]=-np.mean(dfV_shuffled[dfV_shuffled.source == areas[i]].ffb_nc)
        hrt[i]=np.mean(dfV_shuffled[dfV_shuffled.target == areas[i]].ffb_nc)
        hr[i]=(hrs[i]+hrt[i])/2
        hrc[i]=0.5*(-np.mean(dfV_shuffled[dfV_shuffled.source == areas[i]].ffb_c*dfV_shuffled[dfV_shuffled.source == areas[i]].conf)
        +np.mean(dfV_shuffled[dfV_shuffled.target == areas[i]].ffb_c*dfV_shuffled[dfV_shuffled.target == areas[i]].conf))
     
    data=[areas,hrc,hr]
    data=np.transpose(data)
    columns = ['areas','hrc','hr']
    dfi_shuffled = pd.DataFrame(data,columns=columns) 
    
    hr_iter, hrc_iter = iterativeCC(dfi_shuffled,dfV_shuffled,n_iter)
    
    hrc_iter = hrc_iter[['areas',0,n_iter]]
    hr_iter = hr_iter[['areas',0,n_iter]]  
    
    ###########################################################################
    '''Save hierarchy scores of cortical areas in the shuffled data'''
    hr_iter.to_excel(output_dir+'CCshuffled_noconf_iter'+str(i_shuffle)+'.xlsx')
    hrc_iter.to_excel(output_dir+'CCshuffled_conf_iter'+str(i_shuffle)+'.xlsx')
    ###########################################################################
    
    dfV_temp = dfV_shuffled[['source','target','ffb_c','ffb_nc','conf']]
    dfi_temp = hr_iter[['areas',0,n_iter]]
    dfi_temp = dfi_temp.rename(columns={0: "h0", n_iter:"h_iter"})
    dfi_temp_conf = hrc_iter[['areas',0,n_iter]]
    dfi_temp_conf = dfi_temp_conf.rename(columns={0: "h0", n_iter:"h_iter"})

    dfi_t = dfi_temp[['areas','h0','h_iter']]
    dfV_temp = dfV_temp.join(dfi_t.set_index('areas'), on ='source')
    dfV_temp =dfV_temp.rename(columns={"h_iter": "hs_iter"})
    dfV_temp =dfV_temp.rename(columns={"h0": "hs0"})
    dfV_temp = dfV_temp.join(dfi_t.set_index('areas'), on ='target')
    dfV_temp =dfV_temp.rename(columns={"h_iter": "ht_iter"})
    dfV_temp =dfV_temp.rename(columns={"h0": "ht0"})
    dfV_temp = dfV_temp.dropna()
    
    ###########################################################################
    '''global hierarchy score of the shuffled data before & after iteration, without Cre-confidence'''
    hr_init_shuffled[i_shuffle] = np.mean(dfV_temp.ffb_nc*(dfV_temp.ht0 - dfV_temp.hs0))
    hr_iter_shuffled[i_shuffle] = np.mean(dfV_temp.ffb_nc*(dfV_temp.ht_iter - dfV_temp.hs_iter))
    
    dfV_temp = dfV_shuffled[['source','target','ffb_c','ffb_nc','conf']]
    dfi_t = dfi_temp_conf[['areas','h0','h_iter']]
    dfV_temp = dfV_temp.join(dfi_t.set_index('areas'), on ='source')
    dfV_temp =dfV_temp.rename(columns={"h_iter": "hs_iter"})
    dfV_temp =dfV_temp.rename(columns={"h0": "hs0"})
    dfV_temp = dfV_temp.join(dfi_t.set_index('areas'), on ='target')
    dfV_temp =dfV_temp.rename(columns={"h_iter": "ht_iter"})
    dfV_temp =dfV_temp.rename(columns={"h0": "ht0"})
    dfV_temp = dfV_temp.dropna() 
    
     ###########################################################################
    '''global hierarchy score of the shuffled data before & after iteration, with Cre-confidence'''   
    hrc_init_shuffled[i_shuffle] = np.mean(dfV_temp.ffb_c*(dfV_temp.ht0 - dfV_temp.hs0))
    hrc_iter_shuffled[i_shuffle] = np.mean(dfV_temp.ffb_c*(dfV_temp.ht_iter - dfV_temp.hs_iter))
    ###########################################################################    
    
pd.DataFrame(hrc_init_shuffled).to_excel(output_dir +'CC_hg_init_shuffled.xlsx')
pd.DataFrame(hrc_iter_shuffled).to_excel(output_dir +'CC_hg_iter_shuffled.xlsx')


# # Plot global hierarchy scores of 100 shuffled data with the global hierarchy score of the original data


"""Global hierarchy scores of the original cortico-cortical connectivity"""
df_hg_CC = pd.read_excel(input_dir2+'ghs_CC.xlsx')

hg_CC_conf_init = df_hg_CC["hg_CC_conf_init"][0]
hg_CC_conf_iter = df_hg_CC["hg_CC_conf_iter"][0]

hmc_init = (hg_CC_conf_init-np.mean(hrc_init_shuffled))/np.std(hrc_init_shuffled) # Z-score before iteration
hmc_iter = (hg_CC_conf_iter-np.mean(hrc_iter_shuffled))/np.std(hrc_iter_shuffled) # Z-score after iteration

""" Figure showing global hierarchy scores of shuffled data & original CC data before & after iteration """
fig,ax=plt.subplots()
ax.hist(hrc_init_shuffled, bins=10, label='before iterate')
ax.axvline(x=hg_CC_conf_init,linestyle='--')
ax.hist(hrc_iter_shuffled, bins=10, label='after iterate',color='red')
ax.axvline(x=hg_CC_conf_iter,linestyle='--',color='red')
ax.set_xlabel('global hierarchy score',fontsize=16)
ax.set_ylabel('counts',fontsize=16)
ax.set_title('hmc(init)='+str(hmc_init)+'; hmc(iter)='+str(hmc_iter))
ax.legend(loc='upper right')
#fig.savefig(output_dir+"shuffled_globalhierarchy_conf_CC.pdf", bbox_inches='tight')

""" Figure showing global hierarchy scores of shuffled data & original CC data before iteration """
fig,ax=plt.subplots()
ax.hist(hrc_init_shuffled, bins=10, label='before iterate')
ax.axvline(x=hg_CC_conf_init,linestyle='--')
ax.set_xlabel('global hierarchy score',fontsize=16)
ax.set_ylabel('counts',fontsize=16)
ax.set_title('hm(init)='+str(hmc_init))
ax.legend(loc='upper right')
#fig.savefig(output_dir+"shuffled_globalhierarchy_CC_init.pdf", bbox_inches='tight')

""" Figure showing global hierarchy scores of shuffled data & original CC data after iteration """
fig,ax=plt.subplots()
ax.hist(hrc_iter_shuffled, bins=10, label='after iterate',color='red')
ax.axvline(x=hg_CC_conf_iter,linestyle='--',color='red')
ax.set_xlabel('global hierarchy score',fontsize=16)
ax.set_ylabel('counts',fontsize=16)
ax.set_title('hm(iter)='+str(hmc_iter))
ax.legend(loc='upper right')
#fig.savefig(output_dir+"shuffled_globalhierarchy_CC_iter.pdf", bbox_inches='tight')

