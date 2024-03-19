# Retrograde Network Analyses
- Hierarchy searching algorithm based on layer-projection patterns of cortico-cortical, thalamo-cortical, and cortico-thalamic connections constructed from viral tracing experiments using Cre-driver lines
- Code used in "The input-output connectional organization of the cortical-thalamic-basal ganglia circuits and its molecular underpinnings", by Yao et al.
- If you have any questions or suggestions on the code, please contact Hannah Choi (hannahch@gatech.edu) or Eric M Li (eli303@gatech.edu).
- Written 3/18/2024 by Eric M. Li and Hannah Choi.
---
This Python code builds the hierarchy of the mouse cortical and thalamic regions based on clustered layer-projection patterns by maximizing the self-consistency of the obtained hierarchy measured by global hierarchy score, followed by iterative corrections. The hierarchy is first constructed based on anterograde and retrograde cortico-cortical connections only, which can be updated by taking thalamo-cortical and cortico-thalamic connections into account.
## Level of Support
We are releasing this code to the public as a tool we expect others to use. Questions concerning bugs and related issues are welcomed. We expect to address them promptly, pull requests will vetted by our staff before inclusion.

## Hierarchy
### Cortico-cortical hierarchy
#### `run_CC.py`
- The main hierarchy generating code based on anterograde and retrograde cortico-cortical connectivity data.
- Calls:
- - `func_unsupservised_CC_retro.py`
  - `IterativeMethod.py`
- Input files:
- - `AnteroRetro_CC_TC_CT_clusters.xlsx`
- Output files:
- - `CC_conf_iter.xlsx` & `CC_noconf_iter.xlsx` which include hierarchy scores of all cortical regions before & after iterations with or without Cre-confidence. Either of these files is read in `run_TCCT.py` and `func_unsupervised_CT.py`.
  - `inputexpanded_CC.xlsx`, which includes information on CC connections such as sources, targets, and FF/FB mapping. This file is read in `run_TCCT.py` for iteration with TC and CT connections.
  - `ghs_CC.xlsx`, which contains global hierarchy scores of CC connectivity data.
#### `func_unsupervised_CC_retro.py`
- Contains functions needed to map retrograde CC clusters to FF/FB by maximizing the global hierarchy score of the given cortico-cortical connectivity data. Note that anterograde CC clusters have already been optimized in Harris et al (2019) "Hierarchical organization of cortical and thalamic connectivity".
#### `IterativeMethod.py`
- Iterates hierarchy scores of all cortical regions to refine the hierarchy.
- Import `iterativeCC`
### Adding thalamo-cortical and cortico-thalamic hierarchies
The cortico-cortical hierarchy can be updated by including the thalamo-cortical and cortico-thalamic connectivity data. Since the initial thalamic hierarchical positions are determined based on CC hierarchy scores of their source (CT) and target (TC) cortical regions, `run_CC.py` should be called before running `run_TCCT.py`.
#### `run_TCCT.py`
- The main hierarchy generating code based on thalamo-cortical and cortico-thalamic connectivity data.
- Calls:
- - `func_unsupervised_CT.py`
  - `IterativeMethod.py`
- Input files:
- - `AnteroRetro_CC_TC_CT_clusters.xlsx`
  - `CC_conf_iter.xlsx` or `CC_noconf_iter.xlsx`
  - `inputexpanded_CC.xlsx`
- Output files:
- - `TCCT_CCconf_iter.xlsx` (used Cre-confidence for CC hierarchy) or `TCCT_CCnoconf_iter.xlsx` (no Cre-confidence for CC hierarchy), which include hierarchy scores of all cortical and thalamic regions before & after iterations.
  - `inputexpanded_TCCT.xlsx`, which includes information on TC and CT connections such as sources, targets, and FF/FB mappings.
  - `ghs_TCCT.xlsx`, which contains global hierarchy scores of CC+TCCT connectivity data.
#### `func_unsupervised_CT.py`
- Contains functions needed for mapping CT clusters to FF/FB by maximizing the global hierarchy score of the given cortico-thalamic connectivity data. As the initial thalamic hierarchical positions and the FF/FB mapping are based on the cortical hierarchy constructed from `run_CC.py`, the cortical hierarchy scores from CC connectivity should be read accordingly. (In `func_unsupervised_CT.py`, use either `df_cortex=pd.read_excel(r'./Output/CC_conf_iter.xlsx')` or `df_cortex=pd.read_excel(r'./Output/CC_noconf_iter.xlsx')`, depending on whether you use CC hierarchy computed with Cre-confidence or not).
- Note that TC clusters have already been optimized, again in Harris et al (2019) "Hierarchical organization of cortical and thalamic connectivity".
#### `IterativeMethod.py`
- Iterate hierarchy scores of all cortical & thalamic regions
- Import `iterativeTCCT`
#### `IterativeMethod.py`
### Comparison to shuffled connectivity data
To evaluate how hierarchical the mouse brain is, compare the global hierarchy scores of the original connectivity data to the global hierarchy scores of shuffled connectivity data. To do this, use `run_CC_shuffled.py` (with `func_unsupervised_CC_antero.py` and `func_unsupervised_CC_retro.py`) and `run_TCCT_shuffled.py` (with `func_unsupervised_TC_shuffled.py` and `func_unsupervised_CT_shuffled.py`).
### Inter-module and intra-module hierarchy
To find hierarchies of different cortical modules or within each module, run `run_CC_module.py` and `run_TCCT_module.py` to obtain CC and CC+TCCT-based hierarchies. For inter/intra-module hierarchy, instead of searching for the optimal mapping based on the limited numbers of intra-module or inter-module connections every time, use the pre-generated mapping of clusters to either FF or FB direction constructed from all CC, TC, and CT connections, available in `clustermapping.xlsx`. In the paper, we show intra-module hierarchy for the visual module (module = 'Visual') as well as inter-module hierarchy (module ='inter_predefined').

## Centrality
#### `centrality.py`
- Centrality generating code using whole brain CCF-normalized connectivity matrices. Use this code to create correlation plots between centrality and hierarchy scores of brain areas, along with their corresponding Excel files.
- Input files:
- - `combined_ipsi_VN.xlsx`
  - `ipsi_C57_antero_VN.xlsx`
  - `retro_ipsi_VN.xlsx`
- Output files:
- - `antero_inward_centrality.xlsx`, which includes centrality and hierarchy scores for all areas in the hierarchy. Centrality scores are found using the anterograde connectivity matrix with "inward distance" as the measure of distance between two areas.
  - `retro_outward_centrality.xlsx`, which uses the retrograde connectivity matrix and "outward distance" to compute centrality scores of all areas in the hierarchy.
  - `averaged_centrality.xlsx`, which averages the two above centrality scores for each area.
## Modularity
#### `modularity.py`
- Modularity generating code using whole brain connectivity matrices. Use this code to create comparison plots of modularity scores between real and shuffled versions of the matrices, along with modules corresponding to a spatial resolution parameter.
- Input files:
- - `combined_ipsi.xlsx`
  - `antero_ipsi.xlsx`
  - `retro_ipsi.xlsx`
- Output files:
- - `retrograde modularity.svg`, which compares modularity scores between the retrograde connectivity matrix and a shuffled version, over varying values of a spatial resolution parameter (gamma).

## Figures and Results
[@Eric, fix this]The summarized results of pre-generated hierarchy scores and global hierarchy scores are available in the `\results` folder. You can use `figures.py` to produce summary figures based on the pre-generated results.
#### `hierarchy_summary_antero_retro_CreConf.xlsx`
- Contains hierarchy score of all areas at each step of hierarchy generation. Used to generate complete hierarchy plot of cortical areas (`CC iterated`) and all areas (`CC+TCCT iterated`).
- Also contains intra-module and inter-module hierarchy scores, for each step of intra/inter-module hierarchy generation.

#### `gh_comparison_antero_retro_shuffled_CreConf.xlsx`
- Contains global hierarchy scores of shuffled data (x100) and real data for varying steps of hierarchy generation.
- Uses `iter shuffled hg (cortex)`, `iter data hg (cortex)`, `iter shuffled hg (cortex+thal)`, and `iter data hg (cortex+thal)` for global hierarchy score comparison figure. 

## Terms of Use
https://alleninstitute.org/legal/terms-use/
