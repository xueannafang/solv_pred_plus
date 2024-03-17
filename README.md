# solv_pred_plus
 Advanced version of SolvPred. Design solvents by selected features.

log :

- Create plot_corr function to do correlation plot for observable and selected (or all) descriptor. (Mar.15, 2024)
- Create rm_zero function to remove invalid zero data of observables. (Mar.14, 2024)
- Save meta data of load_data part, with full parameters. (Mar.15, 2024)
- Enable hansen distance calculation. (Mar.15, 2024)
- Enable mixture parameters calculation. (Mar.14, 2024)
- Create SolvPredPlus. Complete the load_data function.(Mar.14, 2024)

## How to..

- Inputs:

(1) Experimental data file
```ip_exp_csv_name = "input_exp_data.csv"```

(2) Parameters of all candidate solvents
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

Meta data (plot and csv of explained variance ratio, loading matrix, etc.) will also be saved by default. 

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

All plots will be auto-saved by default. 

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

Full log can also be checked by
```
spp.inv_pca_full_log
```
