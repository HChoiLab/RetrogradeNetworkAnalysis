#!/usr/bin/env python
# coding: utf-8


import logging
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from func_unsupervised_TC_module import fit_TC
from func_unsupervised_CT_module import fit_CT
from IterativeMethod import iterativeTCCT

"""
Created on Sat Mar 2 16:52:07 2024

@author: Eric Li, Hannah Choi
"""
"""
This code finds intra- or inter-module hierarchies of cortical network 
by maximizing the global hierarchy score followed by iterations, 
based on TC and CT connectivity data.
"""


# # Set input and output directories


CreConf = 1                    # 1 if using CC hierarchy with Cre-confidence; 0 if not

input_dir = r'./Input/'        # Directory with the file "AnteroRetro_CC_TC_CT_clusters.xlsx" & "clustermapping.xlsx"
output_dir = r'./Output/module/'      # Directory to save the ouputs from the experimental data


# # Define the module.


module = 'inter_predefined'

# In the paper, we used the following: 'Visual' & 'inter_predefined'
# Possible modules: 'Visual', 'Medial', 'Auditory', 'Somatomotor', 'PFC', 'inter_predefined', 'inter'


# # Clusters of TC and CT source-line-target pairs


xls=pd.ExcelFile(input_dir+"AnteroRetro_CC_TC_CT_clusters.xlsx")

df=pd.read_excel(xls,'anteroCCTC_retroCCCT')

df = df[((df["Target Major Division"] == "isocortex")&(df["Source Major Division"] == "thalamus"))
        |((df["Target Major Division"] == "thalamus")&(df["Source Major Division"] == "isocortex"))] # TC/CT connections


# # Use modules found by retrograde analysis


df.loc[((df["source"] == "GU")|(df["source"] == "VISC")|(df["source"] == "SSs")|(df["source"] == "SSp-bfd")
        |(df["source"] == "SSp-ll")|(df["source"] == "SSp-ul")|(df["source"] == "SSp-un")|(df["source"] == "SSp-n")
        |(df["source"] == "SSp-m")|(df["source"] == "MOp")|(df["source"] == "MOs")),"Cortical Source Module"]='Somatomotor'
df.loc[((df["target"] == "GU")|(df["target"] == "VISC")|(df["target"] == "SSs")|(df["target"] == "SSp-bfd")
        |(df["target"] == "SSp-ll")|(df["target"] == "SSp-ul")|(df["target"] == "SSp-un")|(df["target"] == "SSp-n")
        |(df["target"] == "SSp-m")|(df["target"] == "MOp")|(df["target"] == "MOs")),"Cortical Target Module"]='Somatomotor'

df.loc[((df["source"] == "VISal")|(df["source"] == "VISl")|(df["source"] == "VISp")|(df["source"] == "VISpl")
        |(df["source"] == "VISli")|(df["source"] == "VISpor")|(df["source"] == "VISrl")|(df["source"] == "VISa")
        |(df["source"] == "VISam")|(df["source"] == "VISpm")),"Cortical Source Module"]='Visual'
df.loc[((df["target"] == "VISal")|(df["target"] == "VISl")|(df["target"] == "VISp")|(df["target"] == "VISpl")
        |(df["target"] == "VISli")|(df["target"] == "VISpor")|(df["target"] == "VISrl")|(df["target"] == "VISa")
        |(df["target"] == "VISam")|(df["target"] == "VISpm")),"Cortical Target Module"]='Visual'

df.loc[((df["source"] == "ACAd")|(df["source"] == "ACAv")|(df["source"] == "SSp-tr")|(df["source"] == "RSPagl")
        |(df["source"] == "RSPd")|(df["source"] == "RSPv")),"Cortical Source Module"]='Medial'
df.loc[((df["target"] == "ACAd")|(df["target"] == "ACAv")|(df["target"] == "SSp-tr")|(df["target"] == "RSPagl")
        |(df["target"] == "RSPd")|(df["target"] == "RSPv")),"Cortical Target Module"]='Medial'

df.loc[((df["source"] == "TEa")|(df["source"] == "AUDd")|(df["source"] == "AUDp")|(df["source"] == "AUDpo")
        |(df["source"] == "AUDv")),"Cortical Source Module"]='Auditory'
df.loc[((df["target"] == "TEa")|(df["target"] == "AUDd")|(df["target"] == "AUDp")|(df["target"] == "AUDpo")
        |(df["target"] == "AUDv")),"Cortical Target Module"]='Auditory'

df.loc[((df["source"] == "FRP")|(df["source"] == "PL")|(df["source"] == "ILA")|(df["source"] == "ORBl")
        |(df["source"] == "ORBm")|(df["source"] == "ORBvl")|(df["source"] == "AId")|(df["source"] == "AIv")
        |(df["source"] == "AIp")|(df["source"] == "ECT")|(df["source"] == "PERI")),"Cortical Source Module"]='PFC'
df.loc[((df["target"] == "FRP")|(df["target"] == "PL")|(df["target"] == "ILA")|(df["target"] == "ORBl")
        |(df["target"] == "ORBm")|(df["target"] == "ORBvl")|(df["target"] == "AId")|(df["target"] == "AIv")
        |(df["target"] == "AIp")|(df["target"] == "ECT")|(df["target"] == "PERI")),"Cortical Target Module"]='PFC'

list_module = df["Cortical Target Module"].unique()
list_module = list_module[:list_module.size-1]

clu_ffb=pd.read_excel(input_dir+"clustermapping.xlsx")
clu_ffb = clu_ffb.rename(columns={'TC': 'ffb_TC'})
clu_ffb = clu_ffb.rename(columns={'CT': 'ffb_CT'})


# # For intra-medial, select only the areas within the chosen module


if (module != 'inter') and (module != 'inter_predefined'):
    df = df[(df["Cortical Target Module"] == module)|(df["Cortical Source Module"] == module)]


# # If inter-module, change all the target & source area names to the module name


if (module == 'inter') or (module == 'inter_predefined'): 
    for i_module in range(0,len(list_module)):
        df.loc[df["Cortical Target Module"] == list_module[i_module],'target'] = list_module[i_module]
        df.loc[df["Cortical Source Module"] == list_module[i_module],'source'] = list_module[i_module] 


# # Trim the dataframe


df_antero=df[(df.Antero_Retro == "A")] 
df_antero=df_antero[(df_antero.hemi =="ipsi")&(df_antero.creline != "C57BL/6J / Emx1")&(df_antero.target != "VISC")&(df_antero.source != "SSp-un")&(df_antero.target != "SSp-un")&(df_antero.source != "VISC")] # Consider ipsilateral, Cre connectivity data only
df_antero = df_antero[(df_antero["Target Major Division"] == "isocortex")&(df_antero["Source Major Division"] == "thalamus")] # Consider T-C connections

dfV1 = df_antero[['Antero_Retro','source','target','creline','AnteroCluster']]
dfV1 = dfV1.reset_index(drop=True)
dfV1['AnteroCluster'] = dfV1['AnteroCluster'].apply(np.int64)

df_retro=df[(df.Antero_Retro == "R")] 
df_retro = df_retro[(df_retro["Target Major Division"] == "thalamus")&(df_retro["Source Major Division"] == "isocortex")]  # Consider C-T connections

dfV2 = df_retro[['Antero_Retro','source','target','creline','RetroCluster']]

dfV2 = dfV2.reset_index(drop=True)
dfV2['RetroCluster'] = dfV2['RetroCluster'].apply(np.int64)

frames = [dfV1, dfV2]
dfV = pd.concat(frames)
dfV['AnteroCluster'] = dfV['AnteroCluster'].astype('Int64')
dfV['RetroCluster'] = dfV['RetroCluster'].astype('Int64')

num_antero_clu = len(dfV["AnteroCluster"].unique())-1
num_retro_clu = len(dfV["RetroCluster"].unique())-1

dfVT = dfV.reset_index(drop=True)


# # If inter-module, we may want to find the mapping rule, that is not pre-defined.


if module == 'inter':
    logging.debug("performing initial hierarchy assignment")
    hierarchy_vals_antero = fit_TC(dfV1, module)

    jmax_raw_antero, jmax_antero = np.argmax(hierarchy_vals_antero, axis=0)
    jmax_raw_val_antero = hierarchy_vals_antero[jmax_raw_antero][0]
    jmax_val_antero = hierarchy_vals_antero[jmax_antero][1]
    logging.debug("RESULTS")
    n = num_antero_clu
    logging.debug("(jmax_raw_antero, val) = ({:0{n}b}, {:.3f})".format(jmax_raw_antero, jmax_raw_val_antero, n=n))
    logging.debug("(jmax_antero,     val) = ({:0{n}b}, {:.3f})".format(jmax_antero, jmax_val_antero, n=n))
    
    results = dict(jmax_antero=bin(2**n+jmax_antero),
                   jmax_val_antero=jmax_val_antero,
                   jmax_raw_antero=bin(2**n+jmax_raw_antero),
                   jmax_raw_val_antero=jmax_raw_val_antero)
    hrc_original_antero = jmax_val_antero
    hr_original_antero = jmax_raw_val_antero
    
    print(results)
    hierarchy_vals_retro = fit_CT(dfV2, module)

    jmax_raw_retro, jmax_retro = np.argmax(hierarchy_vals_retro, axis=0)
    jmax_raw_val_retro = hierarchy_vals_retro[jmax_raw_retro][0]
    jmax_val_retro = hierarchy_vals_retro[jmax_retro][1]
    logging.debug("RESULTS")
    n = num_retro_clu
    logging.debug("(jmax_raw_retro, val) = ({:0{n}b}, {:.3f})".format(jmax_raw_retro, jmax_raw_val_retro, n=n))
    logging.debug("(jmax_retro,     val) = ({:0{n}b}, {:.3f})".format(jmax_retro, jmax_val_retro, n=n))
    
    results = dict(jmax_retro=bin(2**n+jmax_retro),
                   jmax_val_retro=jmax_val_retro,
                   jmax_raw_retro=bin(2**n+jmax_raw_retro),
                   jmax_raw_val_retro=jmax_raw_val_retro)
    hrc_original_retro = jmax_val_retro
    hr_original_retro = jmax_raw_val_retro
    
    print(results)


# # Define functions needed.


c0r = 2**(num_retro_clu)
c0a = 2**(num_antero_clu)

def ffb_c_TC (cls):
    """Direction of each TC cluster with confidence"""
    if module == 'inter':
        b=(bin(c0a+jmax_antero)[-(cls)])
        return -(2*int(b)-1) # or -(2*int(b)-1)
    elif module != 'inter':
        return int(clu_ffb[clu_ffb.clu == cls].ffb_TC)   

def ffb_nc_TC (cls):
    """Direction of each TC cluster without confidence"""
    if module == 'inter':
        b=(bin(c0a+jmax_raw_antero)[-(cls)])
        return -(2*int(b)-1)  # or (2*int(b)-1)
    elif module != 'inter':
        return int(clu_ffb[clu_ffb.clu == cls].ffb_TC)
    
def ffb_c_CT (cls):
    """Direction of each CT cluster with confidence"""
    if module == 'inter':
        b=(bin(c0r+jmax_retro)[-(cls)])
        return -(2*int(b)-1) # or -(2*int(b)-1)
    elif module != 'inter':
        return int(clu_ffb[clu_ffb.clu == cls].ffb_CT)   

def ffb_nc_CT (cls):
    """Direction of each CT cluster without confidence"""
    if module == 'inter':
        b=(bin(c0r+jmax_raw_retro)[-(cls)])
        return -(2*int(b)-1)  # or (2*int(b)-1)
    elif module != 'inter':
        return int(clu_ffb[clu_ffb.clu == cls].ffb_CT)

def confidence(df):
    """Returns multiplier which biases towards roughly equal # of FF and FB connections"""    
    count_ff = len(df[df.ffb_c==1])
    count_fb = len(df[df.ffb_c==-1])
    confnew = min(count_ff, count_fb)/(count_ff+count_fb)
    return confnew

def hrf (area):
    '''Hierarchy score of each area without confidence'''
    return (-np.mean(dfVT[dfVT.source == area].ffb_nc)+np.mean(dfVT[dfVT.target == area].ffb_nc))/2

def hrcf (area):
    '''Hierarchy score of each area with confidence'''
    return (-np.mean(dfVT[dfVT.source == area].ffb_c)+np.mean(dfVT[dfVT.target == area].ffb_c))/2


# # Produce expanded data frame with  FF/FB, hierarchy values as source & target for each pair of TC/CT connections


'''ATTENTION! Use the confidence-weighted (biased) mapping of the TC+CT clusters ("ffb_c" & "hrc_s/hrc_t");
Otherwise, the algorithm places all thalamic regions below or above all cortical regions'''

dfVT.loc[(dfVT.Antero_Retro == "R"),"ffb_c"]=dfVT[(dfVT.Antero_Retro == "R")].RetroCluster.apply(ffb_c_CT)
dfVT.loc[(dfVT.Antero_Retro == "R"),"ffb_nc"]=dfVT[(dfVT.Antero_Retro == "R")].RetroCluster.apply(ffb_nc_CT)
dfVT.loc[(dfVT.Antero_Retro == "A"),"ffb_c"]=dfVT[(dfVT.Antero_Retro == "A")].AnteroCluster.apply(ffb_c_TC)
dfVT.loc[(dfVT.Antero_Retro == "A"),"ffb_nc"]=dfVT[(dfVT.Antero_Retro == "A")].AnteroCluster.apply(ffb_nc_TC)

dfVT.loc[:,"hrc_s"]=dfVT["source"].apply(hrcf)
dfVT.loc[:,"hrc_t"]=dfVT["target"].apply(hrcf)
dfVT.loc[:,"hr_s"]=dfVT["source"].apply(hrf)
dfVT.loc[:,"hr_t"]=dfVT["target"].apply(hrf)

conf = confidence(dfVT)

dfVT.to_excel(output_dir+'inputexpanded_TCCT_'+module+'.xlsx')


# # Find hierarchy scores of thalamic areas within a module or of modules


areas = dfV1["source"].unique()
n_areas=len(areas) 
hrs=list(range(0,n_areas))
hrt=list(range(0,n_areas))
hrc=list(range(0,n_areas))

for i in range(0,n_areas):
    hrs[i]=-np.mean(dfVT[dfVT.source == areas[i]].ffb_c)
    if len(dfVT[dfVT.target == areas[i]]) != 0:
        hrt[i]=np.mean(dfVT[dfVT.target == areas[i]].ffb_c)
    else:
        hrt[i]=0
    hrc[i]=(hrs[i]+hrt[i])/2
    
data=[areas,hrc]
data=np.transpose(data)
columns = ['areas','h']
dfiT = pd.DataFrame(data,columns=columns)
dfiT.head()

#dfiT.to_excel(output_dir+'initialhierarchy_TCCT_'+module+'.xlsx')


# # Iterate thalamic + cortical hierarchy scores


n_iter = 20

if CreConf == 1:
    dfiC = pd.read_excel(output_dir+"CC_conf_iter_"+module+".xlsx")
elif CreConf == 0:
    dfiC = pd.read_excel(output_dir+"CC_noconf_iter_"+module+".xlsx")

dfiC['h'] = dfiC[n_iter]
dfVC = pd.read_excel(output_dir+"inputexpanded_CC_"+module+".xlsx")
dfVC = dfVC[["source","target","ffb_nc","ffb_c"]]

if CreConf == 1:
    dfVC["ffb"] = dfVC["ffb_c"]
elif CreConf == 0:
    dfVC["ffb"] = dfVC["ffb_nc"]

dfVT["ffb"] = dfVT["ffb_c"]

dfiT = dfiT[["areas","h"]]
dfiC = dfiC[["areas","h"]]
dfVT = dfVT[["source","target","ffb"]]
dfVC = dfVC[["source","target","ffb"]]

hr_iter = iterativeTCCT(dfiC, dfVC, dfiT, dfVT, n_iter)

iteration=np.arange(0,n_iter+1,1)
n_area=np.shape(hr_iter)[0]
allareas = hr_iter["areas"].unique()

""" Figure of hierarchy score iterations """
fig,ax=plt.subplots()
for i in range(0,n_area):
    y=np.squeeze(np.asarray(hr_iter[hr_iter.areas==allareas[i]].iloc[:,1::1]))
    ax.plot(iteration,y)
ax.set_xlim([0, n_iter])
ax.set_xticks(np.arange(0, n_iter, step=5))
ax.set_xlabel('iter')
ax.set_ylabel('hierarchy value')
ax.set_title('confidence adjusted')
plt.show()

""" Figure showing correlation between hierarchy scores before & after iterations"""
hr_final = hr_iter[:][n_iter]
hr_initial = hr_iter[:][0]   
f = plt.figure()
plt.plot(hr_initial,hr_final,'ro')
plt.xlabel('initial hierarchy (conf)')
plt.ylabel('final hierarchy (conf)')
plt.title('r='+str(np.corrcoef(hr_initial, hr_final)[0, 1]))
plt.show()

'''Save hierarchy scores before and after iteration'''

for i_area in range(0,n_area):
    if hr_iter['areas'][i_area] in list(dfiC['areas']):
        hr_iter.loc[i_area,'CortexThalamus'] = 'C'
    else:
        hr_iter.loc[i_area,'CortexThalamus'] = 'T'

hr_iter = hr_iter[['areas','CortexThalamus', 0,n_iter] ]  
hr_iter_save = hr_iter

if CreConf == 1:
    hr_iter_save.to_excel(output_dir+'TCCT_CCconf_iter_'+module+'.xlsx')
elif CreConf == 0: 
    hr_iter_save.to_excel(output_dir+'TCCT_CCnoconf_iter_'+module+'.xlsx')


#  # Print out global hierarchy scores of CC + TCCT connectivity data before and after iteration


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


print('hg of CC+TC+CT before iterate cortex & thalamus='+str(hg_TCCT_init))
print('hg of CC+TC+CT iterate cortex='+str(hg_cortex_TCCT_iter)) 
print('hg of CC+TC+CT iterate cortex & thalamus='+str(hg_TCCT_iter))

'''Save global hierarchy scores'''
newDF= pd.DataFrame([])
newDF=pd.concat([newDF, pd.DataFrame({'hg_TCCT_init':hg_TCCT_init,
                                 'hg_cortex_TCCT_iter':hg_cortex_TCCT_iter, 'hg_TCCT_iter':hg_TCCT_iter},index=[0])])
newDF.to_excel(output_dir+'ghs_TCCT_'+module+'.xlsx')