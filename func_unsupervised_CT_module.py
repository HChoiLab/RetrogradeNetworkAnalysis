from __future__ import division
from functools import partial
import logging
import numpy as np
import pandas as pd

"""
Created on Sat Mar 2 18:01:52 2024

@author: Eric Li, Hannah Choi
"""
"""
This code contains functions needed for mapping CT clusters to FF/FB by 
maximizing the global hierarchy score of the given cortico-thalamic connectivity data, 
and is called in "run_TCCT_module.py".
"""

def confidence(df):
    """Returns multiplier which prevents all thalamic regions placed below/above cortical regions, by biasing towards equal # of FF and FB connections"""   
    count_ff = len(df[df.ffb==1])
    count_fb = len(df[df.ffb==-1])
    confnew = min(count_ff, count_fb)/(count_ff+count_fb)
 
    return confnew


def hierarchy_values(df, val='ffb'):
    """Returns hierarchy values for thalamic regions as target"""
    target = df.groupby('target')[val].mean()
    levels = target

    return df.target.map(levels) 


def assign_hierarchy_levels(df, j):
    """Returns global hierarchy scores with (hgc)/without (hg) TC confidence 
    (= equal # of FF/FB bias)"""
    def ff_or_fb(cluster):
        b = bin(c0 + j)[-(cluster)]
        return 2*int(b)-1

    n = len(df.RetroCluster.unique())
    c0 = 2**n

    df.loc[:, "ffb"] = df.RetroCluster.apply(ff_or_fb)
    confg = confidence(df)

    ##################################################################
    '''Change the file name accordingly'''
    df_cortex=pd.read_excel(r'./Output/module/CC_conf_iter_inter.xlsx') 
    # df_cortex=pd.read_excel(r'./Output/module/CC_noconf_iter.xlsx') 
    ##################################################################
    
    df = df.join(df_cortex.set_index('areas'), on='source')    
    ##################################################################
    '''Take the last hierarchy scores of cortical regions at the end of iterations (20 steps)'''
    levels_s = df.groupby('source')[20].mean()
    ##################################################################
    hr_t = hierarchy_values(df)
    hrc_t = hr_t
    hr_s = df.source.map(levels_s) 
    hrc_s= hr_s 

    hg = np.mean(df.ffb*(hr_t - hr_s))
    hgc = np.mean(df.ffb*(hrc_t - hrc_s))*confg
    
    logging.debug("{:0{n}b}  (hg, hgc) = ({:.3f}, {:.3f})".format(j, hg, hgc, n=n))

    return hg, hgc


def fit_CT(df, parallel=True, n_procs=-1):
    """iterates through assign_hierarchy_levels with parallel support"""
    def generator():
        n_possible = 2**(len(df.RetroCluster.unique()))
        return range(n_possible)

    hierarchy_vals = []
    func = partial(assign_hierarchy_levels, df)
    if parallel:
        from multiprocessing import Pool, cpu_count
        if n_procs < 1:
            n_procs = cpu_count()
        pool = Pool(processes=n_procs)

        for result in pool.imap(func, generator()):
            hierarchy_vals.append(result)
    else:
        for result in map(func, generator()):
            hierarchy_vals.append(result)

    return hierarchy_vals
