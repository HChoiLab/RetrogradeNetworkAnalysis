#!/usr/bin/env python
# coding: utf-8


import logging
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from func_unsupervised_CT import fit_CT
from IterativeMethod import iterativeTCCT

"""
Published on Sat Mar 2 16:40:37 2024

@author: Eric Li, Hannah Choi
"""
"""
This code finds hierarchy scores of cortical & thalamic areas via iterations, based on
cortico-thalamic (CT) connections and thalamo-cortical (TC) connections, in addition to 
the cortico-cortical (CC) connectivity data.
"""


# # Set input and output directories


input_dir = r'../Input/'     # Directory with the file "AnteroRetro_CC_TC_CT_clusters.xlsx" 
output_dir = r'../Output/'   # Directory to save the ouputs from the experimental data

''' ATTENTION! Change the "df_cortex" accordingly in func_unsupervised_CT as well! '''
CreConf = 1                 # 1 if using CC hierarchy with Cre-confidence; 0 if not


# # Read the excel file with source-target-creline pairs and their cluster numbers. Use TC + CT connections.


xls=pd.ExcelFile(input_dir+"AnteroRetro_CC_TC_CT_clusters.xlsx")
df=pd.read_excel(xls,'anteroCCTC_retroCCCT')

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


# # Map 5 clusters of CT connections to FF/FB directions by maximizing global hierarchy score


hierarchy_vals = fit_CT(dfV2)

jmax_raw_retro, jmax_retro = np.argmax(hierarchy_vals, axis=0)
jmax_raw_val_retro = hierarchy_vals[jmax_raw_retro][0]
jmax_val_retro = hierarchy_vals[jmax_retro][1]
logging.debug("RESULTS")
n = num_retro_clu
logging.debug("(jmax_raw, val) = ({:0{n}b}, {:.3f})".format(jmax_raw_retro, jmax_raw_val_retro, n=n))
logging.debug("(jmax,     val) = ({:0{n}b}, {:.3f})".format(jmax_retro, jmax_val_retro, n=n))

results = dict(jmax_retro=bin(2**n+jmax_retro),
               jmax_val_retro=jmax_val_retro,
               jmax_raw_retro=bin(2**n+jmax_raw_retro),
               jmax_raw_val_retro=jmax_raw_val_retro)
hrc_original = jmax_val_retro         # with Cre-confidence
hr_original = jmax_raw_val_retro      # without Cre-confidence

print(results)


# # Define functions needed


c0r = 2**num_retro_clu

def ff_or_fb_retro_c (cls):
    """Direction of each CT cluster with confidence"""
    b=(bin(c0r+jmax_retro)[-(cls)])
    return -(2*int(b)-1)

def ff_or_fb_retro_nc (cls):
    """Direction of each CT cluster without confidence"""
    b=(bin(c0r+jmax_raw_retro)[-(cls)])
    return (2*int(b)-1) 

def ff_or_fb_antero_c(cls):
    """Direction of each TC cluster with confidence"""
    if cls in {1,3,4,5,7,8,9}:
        return 1
    else:
        return -1

def ff_or_fb_antero_nc(cls):
    """Direction of each TC cluster without confidence"""
    if cls in {1,3,4,5,7,8,9}:
        return 1
    else:
        return -1
    
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


# # Produce expanded data frame with  FF/FB, hierarchy values as source & target for each pair of TC + CT connections


dfVT.loc[(dfVT.Antero_Retro == "R"),"ffb_c"]=dfVT[(dfVT.Antero_Retro == "R")].RetroCluster.apply(ff_or_fb_retro_c)
dfVT.loc[(dfVT.Antero_Retro == "R"),"ffb_nc"]=dfVT[(dfVT.Antero_Retro == "R")].RetroCluster.apply(ff_or_fb_retro_nc)
dfVT.loc[(dfVT.Antero_Retro == "A"),"ffb_c"]=dfVT[(dfVT.Antero_Retro == "A")].AnteroCluster.apply(ff_or_fb_antero_c)
dfVT.loc[(dfVT.Antero_Retro == "A"),"ffb_nc"]=dfVT[(dfVT.Antero_Retro == "A")].AnteroCluster.apply(ff_or_fb_antero_nc)

dfVT.loc[:,"hrc_s"]=dfVT["source"].apply(hrcf)
dfVT.loc[:,"hrc_t"]=dfVT["target"].apply(hrcf)
dfVT.loc[:,"hr_s"]=dfVT["source"].apply(hrf)
dfVT.loc[:,"hr_t"]=dfVT["target"].apply(hrf)

conf = confidence(dfVT)

dfVT.to_excel(output_dir+'inputexpanded_TCCT.xlsx')


# # Find hierarchy score for each of thalamic areas


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


# # Iterate thalamic + cortical hierarchy scores


n_iter = 20

'''Load results from CC hierarchy'''

if CreConf == 0:
    dfiC = pd.read_excel(output_dir+"CC_noconf_iter.xlsx")
elif CreConf == 1:
    dfiC = pd.read_excel(output_dir+"CC_conf_iter.xlsx")
    
dfiC['h'] = dfiC[n_iter]
dfVC = pd.read_excel(output_dir+"inputexpanded_CC.xlsx")
dfVC = dfVC[["source","target","ffb_nc","ffb_c"]]

if CreConf == 0:
    dfVC["ffb"] = dfVC["ffb_nc"]
elif CreConf == 1:
    dfVC["ffb"] = dfVC["ffb_c"]
    
'''Use initial hierarchy scores and cluster directions (FF/FB) 
for TC+CT connections with the "equal #s of FF/FB"-confidence weight'''  

dfVT["ffb"] = dfVT["ffb_c"]

dfiT = dfiT[["areas","h"]]
dfiC = dfiC[["areas","h"]]
dfVT = dfVT[["source","target","ffb"]]
dfVC = dfVC[["source","target","ffb"]]

'''Iterations ''' 
hr_iter = iterativeTCCT(dfiC, dfVC, dfiT, dfVT, n_iter)
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
ax.set_title('confidence adjusted')
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

for i_area in range(0,n_area):
    if hr_iter['areas'][i_area] in list(dfiC['areas']):
        hr_iter.loc[i_area,'CortexThalamus'] = 'C'
    else:
        hr_iter.loc[i_area,'CortexThalamus'] = 'T'

hr_iter = hr_iter[['areas','CortexThalamus', 0,n_iter] ]  

if CreConf == 0:
    hr_iter.to_excel(output_dir+'TCCT_CCnoconf_iter.xlsx')
elif CreConf == 1:
    hr_iter.to_excel(output_dir+'TCCT_CCconf_iter.xlsx')


#  # Print out global hierarchy scores for the TC + CT connectivity data before and after iteration


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
newDF= pd.concat([newDF, pd.DataFrame({'hg_TCCT_init':hg_TCCT_init, 'hg_cortex_TCCT_iter':hg_cortex_TCCT_iter,
                                 'hg_TCCT_iter':hg_TCCT_iter},index=[0])])
newDF.to_excel(output_dir+'ghs_TCCT.xlsx')