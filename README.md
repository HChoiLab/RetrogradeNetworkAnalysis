# Retrograde Network Analyses
- Hierarchy, centrality, and modularity analyses from mouse whole brain connectivity data and matrices.
- Code used in Yao et al (2024) "The input-output connectional organization of the cortical-thalamic-basal ganglia circuits and its molecular underpinnings".
- If you have any questions or suggestions on the code, please contact Hannah Choi (hannahch@gatech.edu) or Eric M. Li (eli303@gatech.edu).
- Updated 10/2/2024 by Eric M. Li. Written 3/19/2024 by Eric M. Li and Hannah Choi.
---
This Python code performs network analyses of the mouse whole brain connectivity constructed based on both anterograde and retrograde Cre-driver line viral tracing experiments. Specifically, this code computes three properties of the mouse whole brain network: hierarchy, centrality, and modularity. The code in the folder `Hierarchy` creates a hierarchy of the mouse cortical, thalamic, and basal ganglia regions by assigning clustered layer-projection patterns to either feedforward or feedback and maximizing a self-consistency global hierarchy score. The hierarchy uses both anterograde and retrograde Cre-driver line viral tracing experiments, beginning with cortico-cortical connections, followed by cortico-thalamic and thalamo-cortical connections, and finished with cortico-basal ganglia connections. 
Using mouse whole brain connectivity matrices, the code `centrality.py` computes closeness centrality scores of the regions. We then create correlation plots between hierarchy and closeness centrality score. 
Finally, the code `modularity.py` creates modules of the mouse cortical brain using the Louvain algorithm by maximizing a "modularity metric", and compares the modularity of the mouse brain to a shuffled version of the mouse connectivity matrix.

## Level of Support
We are making this code available to the public as a helpful tool and guide. We welcome questions concerning bugs and related issues, and expect to address them promptly.

## Hierarchy
### Cortico-cortical hierarchy
#### `run_CC.py`
- The beginning hierarchy code using only anterograde and retrograde cortico-cortical connectivity data.
- Calls:
  - `func_unsupservised_CC_retro.py`
  - `IterativeMethod.py`
- Input files:
  - `AnteroRetro_CC_TC_CT_clusters.xlsx`
- Output files:
  - `CC_conf_iter.xlsx` & `CC_noconf_iter.xlsx` contain hierarchy scores of cortical regions before and after iterations, with or without Cre-confidence. Based on which version is chosen, that Excel file will be read in `run_TCCT.py` and `func_unsupervised_CT.py`.
  - `inputexpanded_CC.xlsx` includes an expanded dataframe on anterograde and retrograde CC connections such as sources, targets, and FF/FB mapping. This file is read in `run_TCCT.py` before iteration with TC and CT connections.
  - `ghs_CC.xlsx` contains overall global hierarchy scores of the hierarchy generated by anterograde and retrograde CC connectivity data.
#### `func_unsupervised_CC_retro.py`
- Import `fit_retro` to map retrograde CC clusters to either FF/FB by maximizing the global hierarchy score of the retrograde CC connectivity data. Note that anterograde CC clusters have already been optimized in Harris et al (2019) "Hierarchical organization of cortical and thalamic connectivity".
#### `IterativeMethod.py`
- Import `iterativeCC` to iterate hierarchy scores of cortical regions after initial creation, refining the hierarchy.
### Adding thalamo-cortical and cortico-thalamic connections to the hierarchy
The CC hierarchy can be updated by including TC (anterograde) and CT (retrograde) connectivity data. Note that initial thalamic hierarchical positions must be based off CC hierarchy scores of their source (CT) and target (TC) cortical regions, so `run_CC.py` should be called prior to running `run_TCCT.py`.
#### `run_TCCT.py`
- Hierarchy code that adds on thalamo-cortical and cortico-thalamic connectivity data.
- Calls:
  - `func_unsupervised_CT.py`
  - `IterativeMethod.py`
- Input files:
  - `AnteroRetro_CC_TC_CT_clusters.xlsx`
  - `CC_conf_iter.xlsx` or `CC_noconf_iter.xlsx`
  - `inputexpanded_CC.xlsx`
- Output files:
  - `TCCT_CCconf_iter.xlsx` (used Cre-confidence for CC hierarchy) or `TCCT_CCnoconf_iter.xlsx` (no Cre-confidence for CC hierarchy) include hierarchy scores of all cortical and thalamic regions before & after iterations.
  - `inputexpanded_TCCT.xlsx` includes an expanded dataframe on TC and CT connections such as sources, targets, and FF/FB mappings.
  - `ghs_TCCT.xlsx` contains overall global hierarchy scores of the entire hierarchy generated by CC+TCCT connectivity data.
#### `func_unsupervised_CT.py`
- Import `fit_CT` to map CT clusters to FF/FB by maximizing the global hierarchy score of the cortico-thalamic connectivity data. Make sure that you use the correct cortical hierarchy scores from CC connectivity data. (In `func_unsupervised_CT.py`, use either `df_cortex=pd.read_excel(r'./Output/CC_conf_iter.xlsx')` or `df_cortex=pd.read_excel(r'./Output/CC_noconf_iter.xlsx')`, depending on whether you use CC hierarchy with Cre-confidence or not).
- Note that TC clusters have already been optimized, again in Harris et al (2019) "Hierarchical organization of cortical and thalamic connectivity".
#### `IterativeMethod.py`
- Import `iterativeTCCT` to iterate hierarchy scores of all cortical & thalamic regions.
### Incorporating connections from cortex to basal ganglia into hierarchy
- Finally, we can incorporate CB (retrograde) connections into the CC+TCCT hierarchy. Note that this is meant to be done after `run_CC.py` and `run_TCCT.py`.
#### `run_CB.py`
- Hierarchy code that incorporates cortico-basal ganglia connectivity data.
- Calls:
  - `IterativeMethod.py`
- Input files:
  - `CB_clusters.csv`
- Output files:
  - `TCCTCBconf_iter.xlsx` (with Cre-confidence for CC hierarchy) or `TCCTCBnoconf_iter.xlsx` (without Cre-confidence for CC hierarchy) include hierarchy scores of all cortical, thalamic, and basal ganglia regions before and after iterations.
  - `inputexpanded_CB.xlsx` contains an expanded dataframe of the CB dataset.
  - `ghs_CB.xlsx` includes global hierarchy scores of the entire CC+TCCT+CB hierarchy.
#### `IterativeMethod.py`
- Imports `iterativeCB` for iterating cortical, thalamic, and basal ganglia areas.

### Comparison to shuffled connectivity data
To determine how hierarchical the mouse brain is, evaluate the lower bound of the global hierarchy score by using a shuffled connectivity dataset. To do this, use `run_CC_shuffled.py` (with `func_unsupervised_CC_antero.py` and `func_unsupervised_CC_retro.py`) and `run_TCCT_shuffled.py` (with `func_unsupervised_TC_shuffled.py` and `func_unsupervised_CT_shuffled.py`). When incorporating thalamic areas into our shuffled hierarchy, there is an additional multiplier that necessitates TC and CT connections have similar FF/FB distributions. Otherwise, we would end up with a hierarchy with all thalamic areas on the top or the bottom, rather than these areas being spread across the hierarchy. Finally, we add CB connections using `run_CB_shuffled.py`. Note that all connections are labeled as feedforward here, consistent with our treatment of the original dataset.
### Inter-module and intra-module hierarchy
To find hierarchies within each cortical module or a hierarchy of cortical modules themselves, run `run_CC_module.py` and `run_TCCT_module.py` to obtain CC and CC+TCCT-based hierarchies. For these hierarchies, instead of searching for the optimal mapping of FF/FB using a very limited dataset of intra-module or inter-module connections, use the pre-generated mapping of clusters found from our hierarchy of all CC, TC, and CT connections, available in `clustermapping.xlsx`. In this paper, we show the inter-module hierarchy (module ='inter_predefined').

## Centrality
#### `centrality.py`
- Centrality analysis code using whole brain CCF-normalized connectivity matrices, based on a Python package NetworkX (https://networkx.org/)
- Use this code to create correlation plots between closeness centrality and hierarchy scores of brain areas, along with their corresponding Excel files.
- Input files:
  - `antero_ipsi_VN.xlsx`
  - `retro_ipsi_VN.xlsx`
- Output files:
  - `antero_inward_centrality.xlsx` includes centrality and hierarchy scores for all areas in the hierarchy. Centrality scores are found using the anterograde connectivity matrix with "inward distance" as the measure of distance between two areas.
  - `retro_outward_centrality.xlsx` uses the retrograde connectivity matrix and "outward distance" to compute centrality scores of all areas in the hierarchy.
  - `averaged_centrality.xlsx` averages the two above centrality scores for each area. This averaged centrality was used for our centrality analysis in this paper.
## Modularity
#### `modularity.py`
- Modularity analysis code using whole mouse brain connectivity matrices, based on Brain Connectivity Toolbox (https://sites.google.com/site/bctnet/)
- Use this code to find cortical modules within the mouse brain at a specified spatial resolution parameter, and create comparison plots of modularity scores between real and shuffled versions of the matrices over this parameter. In the paper, the retrograde connectivity matrix was used. 
- Input files:
  - `combined_ipsi.xlsx`
  - `antero_ipsi.xlsx`
  - `retro_ipsi.xlsx`
## Figures and Results
The summarized results of area hierarchy scores and global hierarchy scores are in the `\results` folder. You can use `figures.py` to produce hierarchy summary figures based on the pre-generated results. Modularity figures comparing real and shuffled connectivity matrices are generated in `modularity.py`, while correlation plots between closeness centrality and hierarchy values can be generated in `centrality.py`.
#### `hierarchy_summary_CreConf.xlsx`
- Contains hierarchy score of all areas at each step of hierarchy generation. Used to generate complete hierarchy plot of the cortical network (`CC iterated`), the cortical-thalamic network (`CC+TCCT iterated`), and the cortical-thalamic-basal ganglia network (`CC+TCCT+CB iterated`).
- Also contains intra-module and inter-module hierarchy scores, for each step of intra/inter-module hierarchy generation.
#### `gh_comparison_shuffled_CreConf.xlsx`
- Contains global hierarchy scores of shuffled data (x100) and real data for varying steps of hierarchy generation.
- Uses `iter shuffled hg (cortex)`, `iter data hg (cortex)`, `iter shuffled hg (cortex+thal)`, `iter data hg (cortex+thal)`, `iter shuffled hg (cortex+thal+bg)`, and `iter data hg (cortex+thal+bg)` for global hierarchy score comparison figure. 

## Terms of Use
https://alleninstitute.org/legal/terms-use/
