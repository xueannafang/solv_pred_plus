# solv_pred_plus
 Advanced version of SolvPred. Design solvents by selected features.

log :

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
```hansen_dist_with = [sm_hsp, prod_hsp]```, 

and claim ```do_hansen_dist``` as ```True``` when calling ```load_exp_data``` function later.



- Create an instance of ```SolvPredPlus```: 

```
spp = SolvPredPlus(all_cand_param_name)
```

- Load experimental data: 

```
spp.load_exp_data(ip_exp_csv_name, do_hansen_dist = True, hansen_dist_with = hansen_dist_with, save_meta_data = True)

```
(The above version will also calculate hansen distance as specified above.)

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

To check further data fed into data processing pipline, call: 

```
spp.non_zero_comb
```