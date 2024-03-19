# SolvPred+
 Advanced version of SolvPred. Inverse dsign of solvents by selected features.

<p>
 <img src=https://github.com/xueannafang/solv_pred_plus/blob/main/sp_plus_2.png width=200>
 </p>

## Background and comparison with previous versions

SolvPred+ is a data-driven inverse solvent design workflow aiming to close the feedback loop for the development of novel materials. 

This project is an extension package for solvent selection toolkit: [SolvPred v2.0](https://github.com/xueannafang/hsp_toolkit_solv_pred_v_2.0) (previous versions: [SolvPred prototype](https://github.com/xueannafang/hsp_toolkit_prototype)). 

In the previous version, solvents are designed by Hansen solubility parameters, only. 

The scope of solvent properties of interest is generalised in this version to fit more complicated situations. 

In this extension version, solvents can be designed by additional features, such as dielectric constant, partition coefficient, etc. 

In addition, first-hand observable information has been incooperated in this workflow, with feauture selection and dimensionalitiy reduction technique (PCA) involved. 

Compared to [SolvPred v2.0](https://github.com/xueannafang/hsp_toolkit_solv_pred_v_2.0), ```SolvPred+``` is object-oriented, with more flexible interface with external input data. 


## How to..

- Load ```SolvPredPlus```:

Start a terminal, make sure ```conda``` environment has been initialised (```conda init```). 

```cd``` to the working directory where ```solv_pred_plus.py``` is stored, then type ```jupyter notebook```.

Either create a new notebook and type

```
%run solv_pred_plus.py
```

, or directly make a copy of the example worksheet ```solv_pred_plus_working_example.ipynb``` ([here](https://github.com/xueannafang/solv_pred_plus/blob/main/solv_pred_plus_working_example.ipynb)) and follow the instruction in the notebook (recommended).


- Inputs:

(1) [Experimental data file](https://github.com/xueannafang/solv_pred_plus/blob/main/input_exp_data.csv)
```ip_exp_csv_name = "input_exp_data.csv"```

(2) [Parameters of all candidate solvents](https://github.com/xueannafang/solv_pred_plus/blob/main/all_candidates_param_at_rt.csv)
```all_cand_param_name "all_candidates_param_at_rt.csv"```

(3) (Optional) Hansen parameters of species to compare with, for example, starting material (s.m.) or product (prod), to calculate Hansen distance with each candidate. The data should follow the following format: 

```
sm_hsp = {
    "comp_name":"sm",
    "D":18.2,
    "P":9.18,
    "H":13.4
}

prod_hsp = {
    "comp_name":"prod",
    "D":17.5,
    "P":5.7,
    "H":2.0
}

```

If Hansen distance needs to be calculated, specify ```hansen_dist_with``` variable: 

```
hansen_dist_with = [sm_hsp, prod_hsp]
```
, 

and claim ```do_hansen_dist``` as ```True``` when calling ```load_exp_data``` function later.



- Create an instance of ```SolvPredPlus``` and load candidate parameter database:

```
spp = SolvPredPlus(all_cand_param_csv, op_folder = "test")
```

If successful, we expect to see

```Loading candidate parameters: success.```

```op_folder``` argeument is the folder name for results, default ```test```.

- Load experimental data: 

```
spp.load_exp_data(ip_exp_csv_name, do_hansen_dist = True, hansen_dist_with = hansen_dist_with, save_meta_data = True)

```
(The above version will also calculate hansen distance as specified above.)

Each entry will return a message looking like: 
```
Process Xf2346 starts..
Calculate hansen distance with sm: success.
Calculate hansen distance with prod: success.
```

If every entry is processed properly, we expect to see

```Loading data: success.```

If meta data has been saved, we will see 

```
op_exp_data_with_full_params.csv saved.
```

[op_exp_data_with_full_params.csv](https://github.com/xueannafang/solv_pred_plus/blob/main/test_1/op_exp_data_with_full_params.csv)

- Clean experimental data and remove invalid group:

```
spp.rm_zero()
```

(In this example, observables equal zero are regarded invalid. User can modify ```rm_zero()``` function if other invalid data or data type needs to be dealt with.)

We can compare the difference before and after cleaning dataset by looking at how many entries are in the working sheet:

```
print(len(spp.all_comb_exp)) #expected as 31 in the example test case
print(len(spp.non_zero_comb)) # expected as 24
```

a successful message should say:
```
Remove invalid/zero data: success.
```


To check further data fed into data processing pipline, call: 

```
spp.non_zero_comb
```
(the dict version)

or 

```
spp.non_zero_comb_df
```
(the dataframe version)


- Generate correlation plot between selected (or all) descriptors with target observable ```Y```: 

```
spp.plot_corr(observable = "Y", plot_all = True, x_to_plot = ["D", "P", "H"], save_plot = True, do_std = False, **kwargs)
```

In ```**kwargs```, we can specify whether to save meta data, by setting ```save_meta_data = True```. 

The correlation matrix of descriptors (only) and descriptors with observables will be saved. 

Note that by default, all the descriptors will be selected and plot out.

Descriptors can alternatively specified in the ```x_to_plot``` argument. (Make sure the descriptor name must match the one fed in the ```all_cand_param_csv``` file and data types are valid.)

We can decide to save figures or not by specifying ```save_plot``` argument. 
Plots will be saved as "corr_plot_" followed by descriptor name in the output results folder.

Example correlation plot: [corr_plot_pi_star.png](https://github.com/xueannafang/solv_pred_plus/blob/main/test_1/corr_plot_pi_star.png) (There will be 10 plots in total; the number of plots equals number of descriptors.)

The ```x``` can be based on original data or standardised data by specifying ```do_std``` argument. 
(Since this part is to look at relationship between individual variables with observable, there's not much difference between std or non-std pattern in the plot.)

The standardisation will be carried out regardless. 

The standardised dataset can be called by:

```
spp.all_std_df
```

where the statistical destails (mean and std) of original dataset fed before standardisation can be checked by:
```
spp.zscore_info
```
The ```zscore_info``` is a dictionary containing numerical column names, mean and std of each column. 


A correlation plot among all descriptors will be presented too. 

Example correlation heatmap for all descriptors: [descriptor_heat_map.png](https://github.com/xueannafang/solv_pred_plus/blob/main/test_1/descriptor_heat_map.png)

The correlation matrix of descriptors (only) can be called by:

```
spp.dscpt_df_corr
```
with zsocred df callable by
```
spp.dscpt_zscored_df
```

The correlation matrix of descriptors and observable can be called by:
```
spp.obs_with_dscpt_df_corr
```
and zsocred df:
```
spp.dscpt_with_obs_zscored_df
```

Example correlation heatmap for all descriptors with observable: [obs_with_descriptor_heat_map.png](https://github.com/xueannafang/solv_pred_plus/blob/main/test_1/obs_with_descriptor_heat_map.png).

Full details of correlation matrix for descriptors: [corr_mat_of_descriptors_only.csv](https://github.com/xueannafang/solv_pred_plus/blob/main/test_1/corr_mat_of_descriptors_only.csv).

Full details of correlation matrix for descriptors and observable: [corr_mat_of_descriptors_with_observable.csv](https://github.com/xueannafang/solv_pred_plus/blob/main/test_1/corr_mat_of_descriptors_with_observable.csv).



- Apply principal component analysis (PCA) and related operations.

To start PCA analysis on full descriptor set:
```
spp.do_pca()
```

The default option will decompose the dataset into equivalent number of principal compoenents (PC), 
i.e., if a ten-dimsensional descriptor set is submitted, 10 PCs will be decomposed. 

In the first round, we recommend to use default settings to evaluate how many PCs to keep (cover ca. 80% explained variance ratio).

If calculation succeeded, we will see
```
PCA starts...
Original descriptors:['D', 'P', 'H', 'epsilon', 'pi_star', 'HBD', 'HBA', 'logP_solv', 'r_sm', 'r_prod']
PCA completed. 
 Please refer to pca_log for details.
```

Full calculation details can be found by:
```
spp.pca_log
```

Meta data (plot and csv of explained variance ratio, loading matrix, etc.) will also be saved by default:

Loading matrix plot: [pca_loading_matrix_n_10.png](https://github.com/xueannafang/solv_pred_plus/blob/main/test_1/pca_loading_matrix_n_10.png).

Full details of loading matrix: [pca_loading_matrix_n_10.csv](https://github.com/xueannafang/solv_pred_plus/blob/main/test_1/pca_loading_matrix_n_10.csv).

Scree plot (explained variance ratio of individual PC): [pca_scree_plot_n_10.png](https://github.com/xueannafang/solv_pred_plus/blob/main/test_1/pca_scree_plot_n_10.png).

Full details of scree plot (explained variance ratio of individual PC): [pca_scree_plot_n_10.csv](https://github.com/xueannafang/solv_pred_plus/blob/main/test_1/pca_scree_plot_n_10.csv).



- Confirm number of PC to keep, generate PCA map

In the test example, first two PCs are sufficient to cover key information of the orginal feature space. 

We can now plot the loading of each feature on the first two PC space. 

(To help with the visualisation, only 2D PCA map will be plot.)

To get projections (loadings) of all descriptors in selected PC space, do

```
spp.get_pc_map(n_comp = 2, with_y = True, save_plot = True, save_meta_data = True)
```

We can modify ```n_comp``` if more than two PCs are required. 

If ```with_y``` is ```True```, the PCA map will be coloured with observable values. 

All plots will be auto-saved by default:

Example PCA map of PC1-PC2 space: [pca_map_with_obs_3d_biplot_PC1_PC2.png](https://github.com/xueannafang/solv_pred_plus/blob/main/test_1/pca_map_with_obs_3d_biplot_PC1_PC2.png).

Descriptor loadings on PC1-PC2 axis: [pca_map_3d_loading_PC1_PC2.png](https://github.com/xueannafang/solv_pred_plus/blob/main/test_1/pca_map_3d_loading_PC1_PC2.png).

Example PCA map of PC2-PC3 space: [pca_map_with_obs_3d_biplot_PC2_PC3.png](https://github.com/xueannafang/solv_pred_plus/blob/main/test_1/pca_map_with_obs_3d_biplot_PC2_PC3.png).

Descriptor loadings on PC2-PC3 axis: [pca_map_3d_loading_PC2_PC3.png](https://github.com/xueannafang/solv_pred_plus/blob/main/test_1/pca_map_3d_loading_PC2_PC3.png).


- Inverse target PC coordinate back to original descriptor space

```
spp.inv_pca_coord(n_comp = 2, pc_to_inv = [5, -1], save_meta_data = True)
```

Make sure ```n_comp``` must match the size of ```pc_to_inv```

If succeeded, descriptors and corresponding value will be presented. 

Meta data will be saved with target coordinate in the file name.

```
inverse_pca_df_5_-1.csv saved.
```

Example inverse pca prediction: [inverse_pca_df_5_-1.csv](https://github.com/xueannafang/solv_pred_plus/blob/main/test_1/inverse_pca_df_5_-1.csv).



Full log can also be checked by
```
spp.inv_pca_full_log
```

(Note that the target PC coordinate, ```pc_to_inv```, in this version, is determined manually. User can select the coordinate based on the trend observed from the PCA map.
This interface will be left open for automatic optimisation, depending on user's requriment.)

- Feature selection by PCA loadings

Now, we can evaluate the contribution of each descriptor to the PC we are interested in..

```
spp.sel_feature(n_comp = 2, pc_to_rank = 'PC1', abs_value = True)
```

```pc_to_rank``` is the PC axis we are looking at.

By default, the absolute contribution, i.e., regardless of symbol, will be compared. (```abs_value = True```)

If ```abs_value``` is set as ```False```, the loading will be compared only based on the exact value calculated in the loading matrix. 

(Note that the strength of contribution is only represented by the absolute value. The symbol only determines postive or negative correlation.)

No matter which mode was executed, the projection of selected PC axis will be stored in ```sel_feature_log```:

```
spp.sel_feature_log
```

The best descriptor will be identified as the one with highest projection on the selected PC axis. 
The descriptor name and correpsonding projection will be presented to user and also stored in the log dictionary. 

For the command above, we expect to see:
```
The absolute projection of all descriptors on PC1 axis is: 

epsilon      0.983794
logP_solv    0.972463
pi_star      0.924214
P            0.864230
D            0.749329
r_prod       0.434760
r_sm         0.408919
HBA          0.266419
HBD          0.175824
H            0.143636
Name: PC1, dtype: float64
The best descriptor is: epsilon. 
 Projection: 0.9837942569191095. 
```


- Revisit the original feature space

Once inverse PCA and feature selection both have been compeleted, we can check what is our feature of interest in the original high dimensional space. 

```
spp.hop_back_to_high_dim(sel_by = 0)
```

By default, the "feature of interest" is determined by the best feature selected in the feature selection step. 

We can alternatively specify the target feature by setting ```sel_by``` as the descriptor name.

For example ```sel_by = "r_sm"```. (Make sure the descriptor name must be valid, i.e., involved in the original dataset.)

The function will return the value of selected feature in the original feature space (i.e., before the PCA transform).

If succeed, we expect to see: 
```
The standardised value of {sel_by} in the original feature space is: {orig_feature_value}
```

For example, 

```
The standardised value of epsilon in the original feature space is: 2.29
```


- Work out coefficients of candidate solvents.

In the final step, we will do the prediction of potential solvent combination that could lead to the target value.

```
spp.pred_mixture(n_solv = 2, target_feature = "default", target_value = "default", pred_mode = "closest", save_meta_data = True, **kwargs)
```

By default, a binary mixture (```n_solv = 2```) based on the best feature selected by PCA (```target_feature``` at ```target_value```) will be predicted. 
The default two candidates are selected by the candidates that are ```closest``` to the ```target_value``` in the candidate set. 
Meta data will be saved automatically:

Two ```.csv``` file containing predicted solvent components and calculation details (including real feature value, deviation from target) will be stored.

Example output for predicted solvent components: [mixture_component_for_epsilon_2.29_2_cand.csv](https://github.com/xueannafang/solv_pred_plus/blob/main/test_1/mixture_component_for_epsilon_2.29_2_cand.csv).

Example output for prediction details (results, error, etc.): [mixture_calculation_info_for_epsilon_2.29_2_cand.csv](https://github.com/xueannafang/solv_pred_plus/blob/main/test_1/mixture_calculation_info_for_epsilon_2.29_2_cand.csv).

If succeeded, we expect to see:

```
['DMSO', 'DMAC'] [2.33598192759192, 1.577976283157882]
First 2 closest candidates: 

solvent
DMSO    0.044486
DMAC    0.713520
Name: epsilon, dtype: float64
['DMSO', 'DMAC'] [2.33598192759192, 1.577976283157882]
Calculating mixture coefficient...

Done.

Coefficient matrix: [array([0.59683364]), array([0.40316636])].

Calculation details: 

{'real_rsl': array([2.03037955]), 'calc_error': array([-0.26111667]), 'percent_error': array([-0.11395029])}
{'cand': ['DMSO', 'DMAC'], 'ratio': [0.5968336404757916, 0.40316635952420843]}
mixture_component_for_epsilon_2.29_2_cand.csv saved.
mixture_calculation_info_for_epsilon_2.29_2_cand.csv saved.
```

The results means a mixture of DMSO (60%) and DMAC (40%) has been constructed. 
The real ```epsilon``` is 2.03, minus 0.26 from the target epsilon: 2.29.

Full results can be found in two csv suggested in the output message. 


An alternative prediction mode is by suggesting solvent candidate manually. 

That can be done by specifying the ```pred_mode``` as ```"specify"```, and feed the solvent name in the kwargs using ```solv_cand = ['solv_name_1', 'solv_name_2']```:

```
spp.pred_mixture(n_solv = 2, pred_mode = "specify", solv_cand = ['DMSO', 'DMAC'])
```

The number of candidates can be modified by varying ```n_solv``` variable. 

Make sure the number suggested in ```solv_cand``` argument must match the number used for ```n_solv```. 

The ```n_solv``` can also be varied in the default mode. 

If ```pred_mode``` is set as ```"full"```, all the possible combinations based on the candidate set will be attempted:

```
spp.pred_mixture(n_solv = 2, target_feature = "default", target_value = "default", pred_mode = "full")
```

If the prediction is successful, we expect to see a list of possible solvent combinations that could led to the target feature value.

The result will be saved with a ```.csv``` file ended by ```all_cand```.
```
mixture_component_for_epsilon_2.29_2_all_cand.csv saved.
```

Example output for prediction based on all possibilties: [mixture_component_for_epsilon_2.29_2_all_cand.csv](https://github.com/xueannafang/solv_pred_plus/blob/main/test_1/mixture_component_for_epsilon_2.29_2_all_cand.csv).

Valid outputs only: [mixture_component_for_epsilon_2.29_2_all_cand_valid.csv](https://github.com/xueannafang/solv_pred_plus/blob/main/test_1/mixture_component_for_epsilon_2.29_2_all_cand_valid.csv)

## All example outputs

If all the ```save_meta_data``` have been specified as ```True```, we expect to see apprx. 26 files in the output folder named by ourselves (when creating the instance of ```SolvPredPlus```):

[test_1](https://github.com/xueannafang/solv_pred_plus/blob/main/test_1)

The exact file number depending on how many descriptors, how many PCs we chose, etc...


## Cite this work

This project is licensed under GPL-3.0.

X. Fang, C. Saunders, E. Gale, N. Fey, C. F.J. Faul, SolvPredPlus - a data-driven inverse solvent design toolkit, 2024, https://github.com/xueannafang/solv_pred_plus.

## Have fun!
