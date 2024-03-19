#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import logging
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from func_unsupervised_CC_retro import fit_retro
from func_unsupervised_CC_antero import fit_antero
from IterativeMethod import iterativeCC

"""
Published on Sat Mar 2 16:48:25 2024

@author: Eric Li, Hannah Choi
"""
"""
This code finds intra- or inter-module hierarchies of cortical network 
by maximizing the global hierarchy score followed by iterations, 
based on the cortico-cortical (CC) connectivity data.
"""


# # Set input and output directories 
# 

# In[ ]:


CreConf = 1 # 1 if using CC hierarchy with Cre-confidence; 0 if not

input_dir = r'./Input/'        # Directory with the file "AnteroRetro_CC_TC_CT_clusters.xlsx" & "clustermapping.xlsx"
output_dir = r'./Output/module/'      # Directory to save the outputs from the experimental data


# # Define the module. 
# 

# In[ ]:


module = 'inter_predefined'

# In the paper, we used the following: 'Visual' & 'inter_predefined'
# Possible modules: 'Visual', 'Medial', 'Auditory', 'Somatomotor', 'PFC', 'inter_predefined', 'inter'


# # Clusters of cortico-cortical source-line-target pairs
# 

# In[ ]:


xls=pd.ExcelFile(input_dir+"AnteroRetro_CC_TC_CT_clusters.xlsx")
df=pd.read_excel(xls,'anteroCCTC_retroCCCT')
df = df[(df["Target Major Division"] == "isocortex")&(df["Source Major Division"] == "isocortex")] # C-C connections


# # Use modules found by retrograde analysis

# In[ ]:


df.loc[((df["source"] == "VISal")|(df["source"] == "VISl")|(df["source"] == "VISp")|(df["source"] == "VISpl")
        |(df["source"] == "VISli")|(df["source"] == "VISpor")|(df["source"] == "VISrl")|(df["source"] == "VISa")
        |(df["source"] == "VISam")|(df["source"] == "VISpm")),"Cortical Source Module"]='Visual'
df.loc[((df["target"] == "VISal")|(df["target"] == "VISl")|(df["target"] == "VISp")|(df["target"] == "VISpl")
        |(df["target"] == "VISli")|(df["target"] == "VISpor")|(df["target"] == "VISrl")|(df["target"] == "VISa")
        |(df["target"] == "VISam")|(df["target"] == "VISpm")),"Cortical Target Module"]='Visual'

df.loc[((df["source"] == "GU")|(df["source"] == "VISC")|(df["source"] == "SSs")|(df["source"] == "SSp-bfd")
        |(df["source"] == "SSp-ll")|(df["source"] == "SSp-ul")|(df["source"] == "SSp-un")|(df["source"] == "SSp-n")
        |(df["source"] == "SSp-m")|(df["source"] == "MOp")|(df["source"] == "MOs")),"Cortical Source Module"]='Somatomotor'
df.loc[((df["target"] == "GU")|(df["target"] == "VISC")|(df["target"] == "SSs")|(df["target"] == "SSp-bfd")
        |(df["target"] == "SSp-ll")|(df["target"] == "SSp-ul")|(df["target"] == "SSp-un")|(df["target"] == "SSp-n")
        |(df["target"] == "SSp-m")|(df["target"] == "MOp")|(df["target"] == "MOs")),"Cortical Target Module"]='Somatomotor'

df.loc[((df["source"] == "FRP")|(df["source"] == "PL")|(df["source"] == "ILA")|(df["source"] == "ORBl")
        |(df["source"] == "ORBm")|(df["source"] == "ORBvl")|(df["source"] == "AId")|(df["source"] == "AIv")
        |(df["source"] == "AIp")|(df["source"] == "ECT")|(df["source"] == "PERI")),"Cortical Source Module"]='PFC'
df.loc[((df["target"] == "FRP")|(df["target"] == "PL")|(df["target"] == "ILA")|(df["target"] == "ORBl")
        |(df["target"] == "ORBm")|(df["target"] == "ORBvl")|(df["target"] == "AId")|(df["target"] == "AIv")
        |(df["target"] == "AIp")|(df["target"] == "ECT")|(df["target"] == "PERI")),"Cortical Target Module"]='PFC'

df.loc[((df["source"] == "ACAd")|(df["source"] == "ACAv")|(df["source"] == "SSp-tr")|(df["source"] == "RSPagl")
        |(df["source"] == "RSPd")|(df["source"] == "RSPv")),"Cortical Source Module"]='Medial'
df.loc[((df["target"] == "ACAd")|(df["target"] == "ACAv")|(df["target"] == "SSp-tr")|(df["target"] == "RSPagl")
        |(df["target"] == "RSPd")|(df["target"] == "RSPv")),"Cortical Target Module"]='Medial'

df.loc[((df["source"] == "TEa")|(df["source"] == "AUDd")|(df["source"] == "AUDp")|(df["source"] == "AUDpo")
        |(df["source"] == "AUDv")),"Cortical Source Module"]='Auditory'
df.loc[((df["target"] == "TEa")|(df["target"] == "AUDd")|(df["target"] == "AUDp")|(df["target"] == "AUDpo")
        |(df["target"] == "AUDv")),"Cortical Target Module"]='Auditory'

list_module = df["Cortical Target Module"].unique()
clu_ffb=pd.read_excel(input_dir+"clustermapping.xlsx")

if CreConf == 1:
    clu_ffb = clu_ffb.rename(columns={'CC_conf_antero': 'ffb_antero'})
    clu_ffb = clu_ffb.rename(columns={'CC_conf_retro': 'ffb_retro'})
elif CreConf == 0:
    clu_ffb = clu_ffb.rename(columns={'CC_noconf_antero': 'ffb_antero'})
    clu_ffb = clu_ffb.rename(columns={'CC_noconf_retro': 'ffb_retro'})


# # For intra-module, select only the areas within the chosen module
# 

# In[ ]:


if (module != 'inter') and (module != 'inter_predefined'):    
    df = df[(df["Cortical Target Module"] == module)&(df["Cortical Source Module"] == module)]


# # If inter-module, change all the target & source area names to the module name
# 

# In[ ]:


if (module == 'inter') or (module == 'inter_predefined'): 
    for i_module in range(0,len(list_module)):
        df.loc[df["Cortical Target Module"] == list_module[i_module],'target'] = list_module[i_module]
        df.loc[df["Cortical Source Module"] == list_module[i_module],'source'] = list_module[i_module] 


# # Trim the dataframe 
# 

# In[ ]:


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

num_clu_antero = len(dfV["AnteroCluster"].unique())-1
num_clu_retro = len(dfV["RetroCluster"].unique())-1

source_areas = dfV["source"].unique()
target_areas = dfV["target"].unique()

dfV_antero = dfV1[["source","target","creline","AnteroCluster", "Antero_Retro"]].copy()
dfV_retro = dfV2[["source","target","creline","RetroCluster", "Antero_Retro"]].copy()


# # If inter-module, we may want to find the mapping rule, that is not pre-defined
# 

# In[ ]:


if module == 'inter':
    logging.debug("performing initial hierarchy assignment")
    hierarchy_vals_antero = fit_antero(dfV1)

    jmax_raw_antero, jmax_antero = np.argmax(hierarchy_vals_antero, axis=0)
    jmax_raw_val_antero = hierarchy_vals_antero[jmax_raw_antero][0]
    jmax_val_antero = hierarchy_vals_antero[jmax_antero][1]
    logging.debug("RESULTS")
    n = len(dfV1.AnteroCluster.unique())
    logging.debug("(jmax_raw_antero, val) = ({:0{n}b}, {:.3f})".format(jmax_raw_antero, jmax_raw_val_antero, n=n))
    logging.debug("(jmax_antero,     val) = ({:0{n}b}, {:.3f})".format(jmax_antero, jmax_val_antero, n=n))
    
    results = dict(jmax_antero=bin(2**n+jmax_antero),
                   jmax_val_antero=jmax_val_antero,
                   jmax_raw_antero=bin(2**n+jmax_raw_antero),
                   jmax_raw_val_antero=jmax_raw_val_antero)
    hrc_original_antero = jmax_val_antero
    hr_original_antero = jmax_raw_val_antero
    print(results)
    
    hierarchy_vals_retro = fit_retro(dfV2)

    jmax_raw_retro, jmax_retro = np.argmax(hierarchy_vals_retro, axis=0)
    jmax_raw_val_retro = hierarchy_vals_retro[jmax_raw_retro][0]
    jmax_val_retro = hierarchy_vals_retro[jmax_retro][1]
    logging.debug("RESULTS")
    n = len(dfV2.RetroCluster.unique())
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
# 

# In[ ]:


c0a=2**num_clu_antero
c0r=2**num_clu_retro

def ffb_c_antero (cls):
    """Direction of each anterograde CC cluster with Cre-confidence"""
    if module == 'inter':
        b=(bin(c0a+jmax_antero)[-(cls)])
        return (2*int(b)-1)
    elif module != 'inter':
        return int(clu_ffb[clu_ffb.clu == cls].ffb_antero)   

def ffb_nc_antero (cls):
    """Direction of each anterograde CC cluster without Cre-confidence"""
    if module == 'inter':
        b=(bin(c0a+jmax_raw_antero)[-(cls)])
        return (2*int(b)-1)
    elif module != 'inter':
        return int(clu_ffb[clu_ffb.clu == cls].ffb_antero)
    
def ffb_c_retro (cls):
    """Direction of each retrograde CC cluster with Cre-confidence"""
    if module == 'inter':
        b=(bin(c0r+jmax_retro)[-(cls)])
        return (2*int(b)-1)
    elif module != 'inter':
        return int(clu_ffb[clu_ffb.clu == cls].ffb_retro)   

def ffb_nc_retro (cls):
    """Direction of each retrograde CC cluster without Cre-confidence"""
    if module == 'inter':
        b=(bin(c0r+jmax_raw_retro)[-(cls)])
        return (2*int(b)-1)
    elif module != 'inter':
        return int(clu_ffb[clu_ffb.clu == cls].ffb_retro)

def cre_confidence1(df):
    """Returns confidence of cre lines"""
    func = lambda x: 1 - np.abs(x.mean())
    return df.groupby('creline')['ffb_c'].transform(func)

def cre_confidence2(df):
    """Returns an alternative confidence of cre lines"""   
    cre_list = df.creline.unique()
    areas = df.source.unique()
    n=len(areas)
    confnew =pd.DataFrame(np.zeros(len(df.index)),index=df.index, columns=['conf'])
 
    for kk in range(0,len(cre_list)):
        df_cre = df[df.creline == cre_list[kk]]
        print(cre_list[kk])
        count_sym = 0
        count_sameffb = 0
        for ii in range(0,n):
            for jj in range(ii+1,n):
                ij_ffb = np.array(df_cre[(df.source == areas[ii])&(df.target == areas[jj])].ffb_c)
                ji_ffb = np.array(df_cre[(df.source == areas[jj])&(df.target == areas[ii])].ffb_c)
                if len(ij_ffb)==1 and len(ji_ffb)==1:
                    count_sym = count_sym+1
                    if ij_ffb == ji_ffb:
                        count_sameffb = count_sameffb+1 
        confnew[df.creline == cre_list[kk]] = 1-count_sameffb/count_sym
    return confnew

def hrf (area):
    '''Hierarchy score of each area or module without Cre-confidence'''
    return (-np.mean(dfV[dfV.source == area].ffb_nc)+np.mean(dfV[dfV.target == area].ffb_nc))/2

def hrcf (area):
    '''Hierarchy score of each area or module with Cre-confidence'''
    return (-np.mean(dfV[dfV.source == area].ffb_c*dfV[dfV.source == area].conf)+np.mean(dfV[dfV.target == area].ffb_c*dfV[dfV.target == area].conf))/2


# # Produce an expanded data frame with Cre-confidence, FF/FB, hierarchy values as source & target for each pair of CC connections
# 

# In[ ]:


dfV.loc[(dfV.Antero_Retro == "R"),"ffb_c"]=dfV[(dfV.Antero_Retro == "R")].RetroCluster.apply(ffb_c_retro)
dfV.loc[(dfV.Antero_Retro == "R"),"ffb_nc"]=dfV[(dfV.Antero_Retro == "R")].RetroCluster.apply(ffb_nc_retro)
dfV.loc[(dfV.Antero_Retro == "A"),"ffb_c"]=dfV[(dfV.Antero_Retro == "A")].AnteroCluster.apply(ffb_c_antero)
dfV.loc[(dfV.Antero_Retro == "A"),"ffb_nc"]=dfV[(dfV.Antero_Retro == "A")].AnteroCluster.apply(ffb_nc_antero)
dfV.loc[:, "conf"] = cre_confidence1(dfV)
dfV.loc[:,"hrc_s"]=dfV["source"].apply(hrcf)
dfV.loc[:,"hrc_t"]=dfV["target"].apply(hrcf)
dfV.loc[:,"hr_s"]=dfV["source"].apply(hrf)
dfV.loc[:,"hr_t"]=dfV["target"].apply(hrf)

dfV.to_excel(output_dir+'inputexpanded_CC_'+module+'.xlsx')


# # Find hierarchy scores of cortical areas within a module or of modules 
# 

# In[ ]:


areas = np.intersect1d(source_areas, target_areas)

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

#dfi.to_excel(output_dir+'initialhierarchy_CC_'+module+'.xlsx')


# # Iterate cortical hierarchy scores to refine the hierarchy levels
# 

# In[ ]:


'''Iterations'''
n_iter = 20
hr_iter, hrc_iter = iterativeCC(dfi,dfV,n_iter)

iteration=np.arange(0,n_iter+1,1)
n_area=len(areas)

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
hrc_iter.to_excel(output_dir+'CC_conf_iter_'+module+'.xlsx')
#hr_iter.to_excel(output_dir+'CC_noconf_iter_'+module+'.xlsx')


# # Print out global hierarchy scores of the CC connectivity data before and after iteration

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
newDF= pd.concat([newDF, pd.DataFrame({'hg_CC_init':hg_CC_init, 'hg_CC_iter':hg_CC_iter,
                                 'hg_CC_conf_init':hg_CC_conf_init, 'hg_CC_conf_iter':hg_CC_conf_iter},index=[0])])
newDF.to_excel(output_dir+'ghs_CC_'+module+'.xlsx')

