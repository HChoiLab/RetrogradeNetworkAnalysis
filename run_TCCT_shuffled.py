#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import logging
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from IterativeMethod import iterativeTCCT
from func_unsupervised_TC_shuffled import fit_TC
from func_unsupervised_CT_shuffled import fit_CT

"""
Published on Sun Mar 10 11:03:14 2024

@author: Eric Li, Hannah Choi
"""
"""
This code generates shuffled versions of the TC+CT connectitivity data and 
computes the global hierarchy scores of the shuffled data.
"""


# # Set input and output directories 
# 

# In[ ]:


input_dir = r'./Input/'        # Directory with the file "AnteroRetro_CC_TC_CT_clusters.xlsx"
input_dir2 = r'./Output/'               # Directory with the file "ghc_TCCT.xls"
output_dir = r'./Output/shuffled/'   # Directory to save the ouputs from the shuffled experimental data

''' ATTENTION! Change the "df_cortex" accordingly in func_unsupervised_TC_shuffled
and func_unsupervised_CT_shuffled as well! '''
CreConf = 1                 # 1 if using CC hierarchy with Cre-confidence; 0 if not


# # Read the excel file with source-target-creline pairs and their cluster numbers. Construct dataframe using only TC+CT connections.
# 

# In[ ]:


xls=pd.ExcelFile(input_dir+"AnteroRetro_CC_TC_CT_clusters.xlsx")
df=pd.read_excel(xls,'anteroCCTC_retroCCCT')

df_antero=df[(df.Antero_Retro == "A")] 
df_antero=df_antero[(df_antero.hemi == "ipsi")&(df_antero.creline != "C57BL/6J / Emx1")&(df_antero.target != "VISC")&(df_antero.source != "SSp-un")&(df_antero.target != "SSp-un")&(df_antero.source != "VISC")]
df_antero = df_antero[(df_antero["Target Major Division"] == "isocortex")&(df_antero["Source Major Division"] == "thalamus")] # Consider T-C connections

dfV1 = df_antero[['Antero_Retro','source','target','creline','AnteroCluster']]
dfV1 = dfV1.reset_index(drop=True)
dfV1['AnteroCluster'] = dfV1['AnteroCluster'].apply(np.int64)

source_areas1 = dfV1["source"].unique()
num_antero_clu = len(dfV1["AnteroCluster"].unique())
total_TC = len(dfV1)

cluster_count_TC = []
for i in range(1, num_antero_clu+1):
    cluster_count_TC.append(np.count_nonzero((dfV1['AnteroCluster'] == i).values))

df_retro=df[(df.Antero_Retro == "R")] 
df_retro = df_retro[(df_retro["Target Major Division"] == "thalamus")&(df_retro["Source Major Division"] == "isocortex")]  # Consider C-T connections

dfV2 = df_retro[['Antero_Retro','source','target','creline','RetroCluster']]
dfV2 = dfV2.reset_index(drop=True)
dfV2['RetroCluster'] = dfV2['RetroCluster'].apply(np.int64)

num_retro_clu = len(dfV2["RetroCluster"].unique())
total_CT = len(dfV2)

cluster_count_CT = []
for i in range(1, num_retro_clu+1):
    cluster_count_CT.append(np.count_nonzero((dfV2['RetroCluster'] == i).values))

dfV_antero = dfV1[["source","target","creline","AnteroCluster", "Antero_Retro"]].copy()
dfV_retro = dfV2[["source","target","creline","RetroCluster", "Antero_Retro"]].copy()


# # Find global hierarchy scores of shuffled TC+CT connectivity data
# 

# In[ ]:


n_iter = 10
n_shuffle = 100
by_creline = 0
line_list_antero = dfV_antero["creline"].unique()
line_list_retro = dfV_retro["creline"].unique()

hr_ct_shuffled = np.zeros(n_shuffle)
hr_iter_c_shuffled = np.zeros(n_shuffle)
hr_iter_ct_shuffled = np.zeros(n_shuffle)

conf_shuffled = np.zeros(n_shuffle)

for i_shuffle in range(0,n_shuffle):
    
    print('i_shuffle='+str(i_shuffle))
    
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
    dfVT_shuffled = pd.concat(frames)
    dfVT_shuffled['AnteroCluster'] = dfVT_shuffled['AnteroCluster'].astype('Int64')
    dfVT_shuffled['RetroCluster'] = dfVT_shuffled['RetroCluster'].astype('Int64')

    source_areas = dfVT_shuffled["source"].unique()
    target_areas = dfVT_shuffled["target"].unique()
    
    areas = source_areas
    n_areas=len(areas)

    dfVT_shuffled = dfVT_shuffled.reset_index(drop=True)
        
    hierarchy_vals_TC = fit_TC(dfV_shuffled_antero,i_shuffle)
    
    jmax_raw_TC, jmax_TC = np.argmax(hierarchy_vals_TC, axis=0)
    jmax_raw_val_TC = hierarchy_vals_TC[jmax_raw_TC][0]
    jmax_val_TC = hierarchy_vals_TC[jmax_TC][1]
    logging.debug("RESULTS")
    n_TC = num_antero_clu
    logging.debug("(jmax_raw_TC, val) = ({:0{n}b}, {:.3f})".format(jmax_raw_TC, jmax_raw_val_TC, n=n_TC))
    logging.debug("(jmax_TC,     val) = ({:0{n}b}, {:.3f})".format(jmax_TC, jmax_val_TC, n=n_TC))

    results_TC = dict(jmax_TC=bin(2**n_TC+jmax_TC),
                   jmax_val_TC=jmax_val_TC,
                   jmax_raw_TC=bin(2**n_TC+jmax_raw_TC),
                   jmax_raw_val_TC=jmax_raw_val_TC)
    
    hierarchy_vals_CT = fit_CT(dfV_shuffled_retro,i_shuffle)
    
    jmax_raw_CT, jmax_CT = np.argmax(hierarchy_vals_CT, axis=0)
    jmax_raw_val_CT = hierarchy_vals_CT[jmax_raw_CT][0]
    jmax_val_CT = hierarchy_vals_CT[jmax_CT][1]
    logging.debug("RESULTS")
    n_CT = num_retro_clu
    logging.debug("(jmax_raw_CT, val) = ({:0{n}b}, {:.3f})".format(jmax_raw_CT, jmax_raw_val_CT, n=n_CT))
    logging.debug("(jmax_CT,     val) = ({:0{n}b}, {:.3f})".format(jmax_CT, jmax_val_CT, n=n_CT))

    results_CT = dict(jmax_CT=bin(2**n_CT+jmax_CT),
                   jmax_val_CT=jmax_val_CT,
                   jmax_raw_CT=bin(2**n_CT+jmax_raw_CT),
                   jmax_raw_val_CT=jmax_raw_val_CT)
    
    multiplier_TC = np.random.choice([-1, 1])
    
    ones_total_TC = 0
    for i in range(1, num_antero_clu+1):    
        ones_total_TC += int(bin(2**num_antero_clu + jmax_TC)[-i]) * cluster_count_TC[i-1]

    ones_total_CT = 0
    for i in range(1, num_retro_clu+1):    
        ones_total_CT += int(bin(2**num_retro_clu + jmax_CT)[-i]) * cluster_count_CT[i-1]

    if ((ones_total_TC > total_TC/2) and (ones_total_CT > total_CT/2)) or ((ones_total_TC < total_TC/2) and (ones_total_CT < total_CT/2)):
        multiplier_CT = multiplier_TC
    else:
        multiplier_CT = -multiplier_TC
        
    c0r=2**n_CT
    c0a=2**n_TC
    
    """Define functions needed"""
    
    def ff_or_fb_retro (cls):
        """Direction of each CT cluster with confidence"""
        b=(bin(c0r+jmax_CT)[-(cls)])
        return multiplier_CT*(2*int(b)-1)
    
    def ff_or_fb_antero(cls):
        """Direction of each TC cluster with confidence"""
        b=(bin(c0a+jmax_TC)[-(cls)])
        return multiplier_TC*(2*int(b)-1)
    
    def ff_or_fb_retro_nc (cls):
        """Direction of each CT cluster without confidence"""
        b=(bin(c0r+jmax_raw_CT)[-(cls)])
        return (2*int(b)-1)
    
    def ff_or_fb_antero_nc(cls):
        """Direction of each TC cluster without confidence"""
        b=(bin(c0a+jmax_raw_TC)[-(cls)])
        return (2*int(b)-1)
    
    def confidence(df):
        """Returns multiplier which biases towards roughly equal # of FF and FB connections"""   
        count_ff = len(df[df.ffb_c==1])
        count_fb = len(df[df.ffb_c==-1])
        confnew = min(count_ff, count_fb)/(count_ff+count_fb)
        return confnew
        
    def hrf (area):
        '''Hierarchy score of each area without confidence'''
        return (-np.mean(dfVT_shuffled[dfVT_shuffled.source == area].ffb_nc)+np.mean(dfVT_shuffled[dfVT_shuffled.target == area].ffb_nc))/2
    
    def hrcf (area):
        '''Hierarchy score of each area with confidence'''
        return (-np.mean(dfVT_shuffled[dfVT_shuffled.source == area].ffb_c)+np.mean(dfVT_shuffled[dfVT_shuffled.target == area].ffb_c))/2

###########################################################################
    
    """Produce an expanded data frame with  FF/FB, hierarchy values as source & target 
    for each pair of TC + CT connections"""
    
    dfVT_shuffled.loc[(dfVT_shuffled.Antero_Retro == "R"),"ffb_c"]=dfVT_shuffled[(dfVT_shuffled.Antero_Retro == "R")].RetroCluster.apply(ff_or_fb_retro)
    dfVT_shuffled.loc[(dfVT_shuffled.Antero_Retro == "R"),"ffb_nc"]=dfVT_shuffled[(dfVT_shuffled.Antero_Retro == "R")].RetroCluster.apply(ff_or_fb_retro_nc)
    dfVT_shuffled.loc[(dfVT_shuffled.Antero_Retro == "A"),"ffb_c"]=dfVT_shuffled[(dfVT_shuffled.Antero_Retro == "A")].AnteroCluster.apply(ff_or_fb_antero)
    dfVT_shuffled.loc[(dfVT_shuffled.Antero_Retro == "A"),"ffb_nc"]=dfVT_shuffled[(dfVT_shuffled.Antero_Retro == "A")].AnteroCluster.apply(ff_or_fb_antero_nc)

    dfVT_shuffled.loc[:,"hrc_s"]=dfVT_shuffled["source"].apply(hrcf)
    dfVT_shuffled.loc[:,"hrc_t"]=dfVT_shuffled["target"].apply(hrcf)
    dfVT_shuffled.loc[:,"hr_s"]=dfVT_shuffled["source"].apply(hrf)
    dfVT_shuffled.loc[:,"hr_t"]=dfVT_shuffled["target"].apply(hrf)

    conf = confidence(dfVT_shuffled)

    ###########################################################################
    '''Finding initial hierarchy score of each thalamic area'''    

    '''Iterate thalamic + cortical hierarchy scores'''
    areas = source_areas1
    n_areas=len(areas) 
    hrs=list(range(0,n_areas))
    hrt=list(range(0,n_areas))
    hrc=list(range(0,n_areas))

    for i in range(0,n_areas):
        hrs[i]=-np.mean(dfVT_shuffled[dfVT_shuffled.source == areas[i]].ffb_c)
        if len(dfVT_shuffled[dfVT_shuffled.target == areas[i]]) != 0:
            hrt[i]=np.mean(dfVT_shuffled[dfVT_shuffled.target == areas[i]].ffb_c)
        else:
            hrt[i]=0
        hrc[i]=(hrs[i]+hrt[i])/2

    data=[areas,hrc]
    data=np.transpose(data)
    columns = ['areas','h']
    dfiT = pd.DataFrame(data,columns=columns)
    ##########################################################
    '''Load results from CC hierarchy'''

    if CreConf == 0:
        dfiC = pd.read_excel(output_dir+'CCshuffled_noconf_iter'+str(i_shuffle)+'.xlsx') 
    elif CreConf == 1:
        dfiC = pd.read_excel(output_dir+'CCshuffled_conf_iter'+str(i_shuffle)+'.xlsx')  

    dfiC['h'] = dfiC[10]
    dfVC = pd.read_excel(output_dir+'inputexpanded_CC_shuffled'+str(i_shuffle)+'.xlsx')
    dfVC = dfVC[["source","target","ffb_nc","ffb_c"]]

    if CreConf == 0:
        dfVC["ffb"] = dfVC["ffb_nc"]
    elif CreConf == 1:
        dfVC["ffb"] = dfVC["ffb_c"]
    ##########################################################    
    
    
    dfVT_shuffled["ffb"] = dfVT_shuffled["ffb_c"]
    dfiT = dfiT[["areas","h"]]
    dfiC = dfiC[["areas","h"]]
    dfVT = dfVT_shuffled[["source","target","ffb"]]
    dfVC = dfVC[["source","target","ffb"]]
    
    hr_iter = iterativeTCCT(dfiC, dfVC, dfiT, dfVT, n_iter)    
    
    iteration=np.arange(0,n_iter+1,1)
    n_area=np.shape(hr_iter)[0]
    allareas = hr_iter["areas"].unique()
    
    for i_area in range(0,n_area):
        if hr_iter['areas'][i_area] in list(dfiC['areas']):
            hr_iter.loc[i_area,'CortexThalamus'] = 'C'
        else:
            hr_iter.loc[i_area,'CortexThalamus'] = 'T'
    
    hr_iter = hr_iter[['areas','CortexThalamus', 0,n_iter] ]  
   
    if CreConf == 1:
        hr_iter.to_excel(output_dir+'TCCT_CCconf_iter'+str(i_shuffle)+'.xlsx')
    elif CreConf == 0:
        hr_iter.to_excel(output_dir+'TCCT_CCnoconf_iter'+str(i_shuffle)+'.xlsx')
    ###########################################################################


    ###########################################################################
    '''global hierarchy score of the shuffled data before & after iteration'''
    
    dfi_TCCT = hr_iter[["CortexThalamus","areas",0,n_iter]]
    dfi_TCCT = dfi_TCCT.rename(columns={0: "h0", n_iter:"h_iter"})
    
    dfV_CC = dfVC[['source','target','ffb']]
    dfV_TCCT = dfVT[["source","target","ffb"]]
    
    
    dfi_cortex1 = dfi_TCCT[(dfi_TCCT.CortexThalamus == 'C')]
    dfi_cortex1 = dfi_cortex1[['areas','h_iter']]
    dfV_CC = dfV_CC.join(dfi_cortex1.set_index('areas'), on ='source')
    dfV_CC=dfV_CC.rename(columns={"h_iter": "hs"})
    dfV_CC = dfV_CC.join(dfi_cortex1.set_index('areas'), on ='target')
    dfV_CC=dfV_CC.rename(columns={"h_iter": "ht"})
    dfV_CC = dfV_CC.dropna() 
    hg_CC_1 = dfV_CC.ffb*(dfV_CC.ht- dfV_CC.hs)
    
    dfi_thalamus1=dfi_TCCT[(dfi_TCCT.CortexThalamus == 'T')]
    dfi_thalamus1 = dfi_thalamus1[['areas','h_iter']]
    dfV_TCCT = dfV_TCCT.join(dfi_thalamus1.set_index('areas'), on ='source')
    dfV_TCCT=dfV_TCCT.rename(columns={"h_iter": "hs"})
    dfV_TCCT = dfV_TCCT.join(dfi_cortex1.set_index('areas'), on ='target')
    dfV_TCCT=dfV_TCCT.rename(columns={"h_iter": "ht"})
    dfV_TCCT = dfV_TCCT.dropna() 
    hg_TCCT_1 = dfV_TCCT.ffb*(dfV_TCCT.ht- dfV_TCCT.hs)
    
    hg_cortex_TCCT_iter = np.mean(hg_CC_1)
    hg_TCCT_iter = np.mean(pd.concat([hg_CC_1, hg_TCCT_1]))
    
    dfV_CC = dfVC[['source','target','ffb']]
    dfV_TCCT = dfVT[["source","target","ffb"]]
    
    dfi_cortex1 = dfi_TCCT[(dfi_TCCT.CortexThalamus == 'C')]
    dfi_cortex1 = dfi_cortex1[['areas','h0']]
    dfV_CC = dfV_CC.join(dfi_cortex1.set_index('areas'), on ='source')
    dfV_CC=dfV_CC.rename(columns={"h0": "hs"})
    dfV_CC = dfV_CC.join(dfi_cortex1.set_index('areas'), on ='target')
    dfV_CC=dfV_CC.rename(columns={"h0": "ht"})
    dfV_CC = dfV_CC.dropna() 
    hg_CC_1 = dfV_CC.ffb*(dfV_CC.ht- dfV_CC.hs)
    
    dfi_thalamus1=dfi_TCCT[(dfi_TCCT.CortexThalamus == 'T')]
    dfi_thalamus1 = dfi_thalamus1[['areas','h0']]
    dfV_TCCT = dfV_TCCT.join(dfi_thalamus1.set_index('areas'), on ='source')
    dfV_TCCT=dfV_TCCT.rename(columns={"h0": "hs"})
    dfV_TCCT = dfV_TCCT.join(dfi_cortex1.set_index('areas'), on ='target')
    dfV_TCCT=dfV_TCCT.rename(columns={"h0": "ht"})
    dfV_TCCT = dfV_TCCT.dropna() 
    hg_TCCT_1 = dfV_TCCT.ffb*(dfV_TCCT.ht- dfV_TCCT.hs)
    
    hg_cortex_TCCT_init = np.mean(hg_CC_1)
    hg_TCCT_init  = np.mean(pd.concat([hg_CC_1, hg_TCCT_1]))
    
  
    hr_ct_shuffled[i_shuffle] = hg_TCCT_init
    hr_iter_c_shuffled[i_shuffle] = hg_cortex_TCCT_iter
    hr_iter_ct_shuffled[i_shuffle] = hg_TCCT_iter

pd.DataFrame(hr_ct_shuffled).to_excel(output_dir+'shuffled_hg_TCCT_init_cortexthalamus.xlsx')
pd.DataFrame(hr_iter_c_shuffled).to_excel(output_dir+'shuffled_hg_TCCT_iter_cortex.xlsx')
pd.DataFrame(hr_iter_ct_shuffled).to_excel(output_dir+'shuffled_hg_TCCT_iter_cortexthalamus.xlsx')


# # Plot global hierarchy scores of 100 shuffled data with the global hierarchy score of the original data
# 

# In[ ]:


"""Global hierarchy scores of the original cortico-thalamic + thalamo-cortical connectivity"""

df_hg_TC = pd.read_excel(input_dir2+'ghs_TCCT.xlsx')

hg_all_init = df_hg_TC["hg_TCCT_init"][0]
hg_cortex_iter = df_hg_TC["hg_cortex_TCCT_iter"][0]
hg_all_iter = df_hg_TC["hg_TCCT_iter"][0] 

hm1 = (hg_all_init-np.mean(hr_ct_shuffled))/np.std(hr_ct_shuffled)              # Z-score for thalamus+cortex before iteration
hm2 = (hg_cortex_iter-np.mean(hr_iter_c_shuffled))/np.std(hr_iter_c_shuffled)   # Z-score for cortex after iteration
hm3 = (hg_all_iter-np.mean(hr_iter_ct_shuffled))/np.std(hr_iter_ct_shuffled)    # Z-score for thalamus+cortex after iteration

""" Figure showing global hierarchy scores of shuffled data & original CT+TC data before iteration, for cortex+thalamus """
fig,ax=plt.subplots()
bins=25
ax.hist(hr_ct_shuffled, bins=bins, label='confidence adjusted')
ax.axvline(x=hg_all_init,linestyle='--')
ax.set_xlabel('global hierarchy score',fontsize=16)
ax.set_ylabel('counts',fontsize=16)
ax.set_title('hm='+str(hm1))
ax.legend(loc='upper right')
#fig.savefig(output_dir+"shuffledgh_TCCT_all_init.pdf", bbox_inches='tight')

""" Figure showing global hierarchy scores of shuffled data & original CT+TC data after iteration, for cortex only """
fig,ax=plt.subplots()
bins=25
ax.hist(hr_iter_c_shuffled, bins=bins, label='confidence adjusted')
ax.axvline(x=hg_cortex_iter,linestyle='--')
ax.set_xlabel('global hierarchy score',fontsize=16)
ax.set_ylabel('counts',fontsize=16)
ax.set_title('hm='+str(hm2))
ax.legend(loc='upper right')
#fig.savefig(output_dir+"shuffledgh_TCCT_cortex_iter.pdf", bbox_inches='tight')

""" Figure showing global hierarchy scores of shuffled data & original CT+TC data after iteration, for cortex+thalamus """
fig,ax=plt.subplots()
bins=25
ax.hist(hr_iter_ct_shuffled, bins=bins, label='confidence adjusted')
ax.axvline(x=hg_all_iter,linestyle='--')
ax.set_xlabel('global hierarchy score',fontsize=16)
ax.set_ylabel('counts',fontsize=16)
ax.set_title('hm='+str(hm3))
ax.legend(loc='upper right')
#fig.savefig(output_dir+"shuffledgh_TCCT_all_iter.pdf", bbox_inches='tight')

