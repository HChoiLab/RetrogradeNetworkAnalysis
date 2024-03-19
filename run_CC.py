#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import logging
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from func_unsupervised_CC_retro import fit_retro
from IterativeMethod import iterativeCC

"""
Published on Sat Mar 2 16:20:24 2024

@author: Eric Li, Hannah Choi
"""
"""
This code finds hierarchy scores of cortical areas by maximizing the global hierarchy score 
followed by iterations, based on the cortico-cortical (CC) connectivity data from anterograde
and retrograde mouse brain connections.
"""


# # Set input and output directories

# In[ ]:


input_dir = r'./Input/'        # Directory with the file "AnteroRetro_CC_TC_CT_clusters.xlsx"
output_dir = r'./Output/'      # Directory to save the ouputs from the experimental data


# # Read the excel file with source-target-creline pairs and their cluster numbers. Use only the CC connections.

# In[ ]:


xls=pd.ExcelFile(input_dir+"AnteroRetro_CC_TC_CT_clusters.xlsx")

df=pd.read_excel(xls,'anteroCCTC_retroCCCT')

df_antero=df[(df.Antero_Retro == "A")] 
df_antero=df_antero[(df_antero.hemi =="ipsi")&(df_antero.creline != "C57BL/6J / Emx1")&(df_antero.target != "VISC")&(df_antero.source != "SSp-un")&(df_antero.target != "SSp-un")&(df_antero.source != "VISC")] # Consider ipsilateral, Cre connectivity data only
df_antero = df_antero[(df_antero["Target Major Division"] == "isocortex")&(df_antero["Source Major Division"] == "isocortex")] # Consider C-C connections only

dfV1 = df_antero[['Antero_Retro','source','target','creline','AnteroCluster']]
dfV1 = dfV1.reset_index(drop=True)
dfV1['AnteroCluster'] = dfV1['AnteroCluster'].apply(np.int64)

df_retro=df[(df.Antero_Retro == "R")] 
df_retro = df_retro[(df_retro["Target Major Division"] == "isocortex")&(df_retro["Source Major Division"] == "isocortex")]  # Consider C-C connections only

dfV2 = df_retro[['Antero_Retro','source','target','creline','RetroCluster']]
dfV2 = dfV2.reset_index(drop=True)
dfV2['RetroCluster'] = dfV2['RetroCluster'].apply(np.int64) 

frames = [dfV1, dfV2]
dfV = pd.concat(frames)
dfV['AnteroCluster'] = dfV['AnteroCluster'].astype('Int64')
dfV['RetroCluster'] = dfV['RetroCluster'].astype('Int64')

num_antero_clu = len(dfV["AnteroCluster"].unique())-1
num_retro_clu = len(dfV["RetroCluster"].unique())-1

dfV = dfV.reset_index(drop=True)


# # Map 11 clusters of retrograde CC connections to FF/FB directions by maximizing the global hierarchy score

# In[ ]:


hierarchy_vals = fit_retro(dfV2)

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

# In[ ]:


c0r = 2**num_retro_clu

def ff_or_fb_retro_c (cls):
    """Direction of each retrograde CC cluster with Cre-confidence"""
    b=(bin(c0r+jmax_retro)[-(cls)])
    return (2*int(b)-1)    # or -(2*int(b)-1), depending on the final direction of the hierarchy

def ff_or_fb_retro_nc (cls):
    """Direction of each retrograde CC cluster without Cre-confidence"""
    b=(bin(c0r+jmax_raw_retro)[-(cls)])
    return (2*int(b)-1)

def ff_or_fb_antero_c(cls):
    """Direction of each anterograde CC cluster with Cre-confidence"""
    if cls in {1,3,4,5,7,8}:
        return 1
    else:
        return -1

def ff_or_fb_antero_nc(cls):
    """Direction of each anterograde CC cluster without Cre-confidence"""
    if cls in {1,3,4,5,7,8}:
        return 1
    else:
        return -1

def cre_confidence1(df):
    """Returns confidence of Cre lines"""
    func = lambda x: 1 - np.abs(x.mean())
    return df.groupby('creline')['ffb_c'].transform(func)

def hrf (area):
    '''Hierarchy score of each area without Cre-confidence'''
    return (-np.mean(dfV[dfV.source == area].ffb_nc)+np.mean(dfV[dfV.target == area].ffb_nc))/2

def hrcf (area):
    '''Hierarchy score of each area with Cre-confidence'''
    return (-np.mean(dfV[dfV.source == area].ffb_c*dfV[dfV.source == area].conf)+np.mean(dfV[dfV.target == area].ffb_c*dfV[dfV.target == area].conf))/2


# # Produce expanded data frame with FF/FB hierarchy values as source & target for each pair of CC connections

# In[ ]:


dfV.loc[(dfV.Antero_Retro == "R"),"ffb_c"]=dfV[(dfV.Antero_Retro == "R")].RetroCluster.apply(ff_or_fb_retro_c)
dfV.loc[(dfV.Antero_Retro == "R"),"ffb_nc"]=dfV[(dfV.Antero_Retro == "R")].RetroCluster.apply(ff_or_fb_retro_nc)
dfV.loc[(dfV.Antero_Retro == "A"),"ffb_c"]=dfV[(dfV.Antero_Retro == "A")].AnteroCluster.apply(ff_or_fb_antero_c)
dfV.loc[(dfV.Antero_Retro == "A"),"ffb_nc"]=dfV[(dfV.Antero_Retro == "A")].AnteroCluster.apply(ff_or_fb_antero_nc)
dfV.loc[:, "conf"] = cre_confidence1(dfV)
dfV.loc[:,"hrc_s"]=dfV["source"].apply(hrcf)
dfV.loc[:,"hrc_t"]=dfV["target"].apply(hrcf)
dfV.loc[:,"hr_s"]=dfV["source"].apply(hrf)
dfV.loc[:,"hr_t"]=dfV["target"].apply(hrf)

dfV.to_excel(output_dir+'inputexpanded_CC.xlsx')


# # Find hierarchy score for each of cortical areas

# In[ ]:


source_areas = dfV1["source"].unique()
target_areas = dfV2["target"].unique()
areas = np.concatenate ([source_areas,target_areas],axis=None)
areas = np.unique(areas)
n_areas=len(areas)
hrs=list(range(0,n_areas))
hrt=list(range(0,n_areas))
hr=list(range(0,n_areas))
hrc=list(range(0,n_areas))

for i in range(0,n_areas):
    hrs[i]=-np.mean(dfV[dfV.source == areas[i]].ffb_nc)
    hrt[i]=np.mean(dfV[dfV.target == areas[i]].ffb_nc)
    hr[i]=(hrs[i]+hrt[i])/2
    hrc[i]=0.5*(-np.mean(dfV[dfV.source == areas[i]].ffb_c*dfV[dfV.source == areas[i]].conf)+np.mean(dfV[dfV.target == areas[i]].ffb_c*dfV[dfV.target == areas[i]].conf))
    
data=[areas,hrc,hr]
data=np.transpose(data)
columns = ['areas','hrc','hr']
dfi = pd.DataFrame(data,columns=columns)
dfi.head()

dfi.to_excel(output_dir+'initialhierarchy_CC.xlsx')


# # Iterate cortical hierarchy scores

# In[ ]:


""" Iterations """

n_iter = 20
hr_iter, hrc_iter = iterativeCC(dfi,dfV,n_iter)

iteration=np.arange(0,n_iter+1,1)
n_area=len(source_areas)

""" Figure of hierarchy score iterations with Cre-confidence """

fig,ax=plt.subplots()
for i in range(0,n_area):
    y=np.squeeze(np.asarray(hrc_iter[hrc_iter.areas==areas[i]].iloc[:,1::1]))
    ax.plot(iteration,y)

ax.set_xlim([0, n_iter])
ax.set_xticks(np.arange(0, n_iter, step=5))
ax.set_xlabel('iter')
ax.set_ylabel('hierarchy value')
ax.set_title('confidence adjusted')
plt.show()

""" Figure of hierarchy score iterations without Cre-confidence """

fig,ax=plt.subplots()
for i in range(0,n_area):
    y=np.squeeze(np.asarray(hr_iter[hr_iter.areas==areas[i]].iloc[:,1::1]))
    ax.plot(iteration,y)

ax.set_xlim([0, n_iter])
ax.set_xticks(np.arange(0, n_iter, step=5))
ax.set_xlabel('iter')
ax.set_ylabel('hierarchy value')
ax.set_title('confidence not adjusted')
plt.show()

""" Figure showing correlation between hierarchy scores before & after iterations with Cre-confidence """

hrc_final = hrc_iter[:][n_iter]
hrc_initial = hrc_iter[:][0]    
f = plt.figure()
plt.plot(hrc_initial,hrc_final,'ro')
plt.xlabel('initial hierarchy (conf)')
plt.ylabel('final hierarchy (conf)')
plt.title('r='+str(np.corrcoef(hrc_initial, hrc_final)[0, 1]))
plt.show()  

""" Figure showing correlation between hierarchy scores before & after iterations without Cre-confidence """

hr_final = hr_iter[:][n_iter]
hr_initial = hr_iter[:][0]   
f = plt.figure()
plt.plot(hr_initial,hr_final,'ro')
plt.xlabel('initial hierarchy (noconf)')
plt.ylabel('final hierarchy (noconf)')
plt.title('r='+str(np.corrcoef(hr_initial, hr_final)[0, 1]))
plt.show()

'''Save hierarchy scores before and after iteration'''
hrc_iter = hrc_iter[['areas',0,n_iter]]
hr_iter = hr_iter[['areas',0,n_iter] ]  
hrc_iter.to_excel(output_dir+'CC_conf_iter.xlsx')  # Save before & after iteration hierarchy scores with Cre-confidence
hr_iter.to_excel(output_dir+'CC_noconf_iter.xlsx') # Save before & after iteration hierarchy scores without Cre-confidence


# # Print out global hierarchy scores for the CC connectivity data before and after iteration

# In[ ]:


dfV_temp = dfV[['source','target','ffb_c','ffb_nc','conf']]
dfi_temp = hr_iter[['areas',0,n_iter]]
dfi_temp_conf = hrc_iter[['areas',0,n_iter]]
dfi_temp = dfi_temp.rename(columns={0: "h0", n_iter:"h_iter"})
dfi_temp_conf = dfi_temp_conf.rename(columns={0: "h0", n_iter:"h_iter"})

dfi_t = dfi_temp[['areas','h0']]
dfV_temp = dfV_temp.join(dfi_t.set_index('areas'), on ='source')
dfV_temp =dfV_temp.rename(columns={"h0": "hs_init"})
dfV_temp = dfV_temp.join(dfi_t.set_index('areas'), on ='target')
dfV_temp =dfV_temp.rename(columns={"h0": "ht_init"})
dfV_temp = dfV_temp.dropna()

dfi_t = dfi_temp[['areas','h_iter']]
dfV_temp = dfV_temp.join(dfi_t.set_index('areas'), on ='source')
dfV_temp =dfV_temp.rename(columns={"h_iter": "hs_iter"})
dfV_temp = dfV_temp.join(dfi_t.set_index('areas'), on ='target')
dfV_temp =dfV_temp.rename(columns={"h_iter": "ht_iter"})
dfV_temp = dfV_temp.dropna()

hg_CC_init = np.mean(dfV_temp.ffb_nc*(dfV_temp.ht_init - dfV_temp.hs_init))
hg_CC_iter = np.mean(dfV_temp.ffb_nc*(dfV_temp.ht_iter - dfV_temp.hs_iter))


dfV_temp = dfV[['source','target','ffb_c','ffb_nc','conf']]
dfi_t= dfi_temp_conf[['areas','h0']]
dfV_temp= dfV_temp.join(dfi_t.set_index('areas'), on ='source')
dfV_temp =dfV_temp.rename(columns={"h0": "hs_init"})
dfV_temp = dfV_temp.join(dfi_t.set_index('areas'), on ='target')
dfV_temp =dfV_temp.rename(columns={"h0": "ht_init"})
dfV_temp = dfV_temp.dropna()

dfi_t = dfi_temp_conf[['areas','h_iter']]
dfV_temp = dfV_temp.join(dfi_t.set_index('areas'), on ='source')
dfV_temp =dfV_temp.rename(columns={"h_iter": "hs_iter"})
dfV_temp = dfV_temp.join(dfi_t.set_index('areas'), on ='target')
dfV_temp =dfV_temp.rename(columns={"h_iter": "ht_iter"})
dfV_temp = dfV_temp.dropna()

hg_CC_conf_init = np.mean(dfV_temp.ffb_c*(dfV_temp.ht_init - dfV_temp.hs_init))
hg_CC_conf_iter = np.mean(dfV_temp.ffb_c*(dfV_temp.ht_iter - dfV_temp.hs_iter))

print('hg of CC cortex='+str(hg_CC_init))
print('hg of CC iterate cortex='+str(hg_CC_iter))

print('hgc of CC cortex='+str(hg_CC_conf_init))
print('hgc of CC iterate cortex='+str(hg_CC_conf_iter))

'''Save global hierarchy scores'''
newDF= pd.DataFrame([])
newDF=pd.concat([newDF, pd.DataFrame({'hg_CC_init':hg_CC_init, 'hg_CC_iter':hg_CC_iter,
                                 'hg_CC_conf_init':hg_CC_conf_init, 'hg_CC_conf_iter':hg_CC_conf_iter},index=[0])])
newDF.to_excel(output_dir+'ghs_CC.xlsx')

