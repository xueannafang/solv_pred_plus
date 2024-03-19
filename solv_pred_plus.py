import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import zscore
from matplotlib.colors import ListedColormap
from sklearn.decomposition import PCA
import os
from itertools import combinations


class SolvPredPlus:
    
    def __init__(self, all_cand_param_csv, op_folder = "test", **kwargs):
        self.load_cand_param(all_cand_param_csv, **kwargs)
        
        if os.path.isdir(op_folder):
            pass
        else:
            os.mkdir(op_folder)
        self.op_folder_path = os.getcwd() + ".\\" + op_folder
        
    
    def load_cand_param(self, all_cand_param_csv, **kwargs):
        """
        load database for all solvent parameters for candidates
        all_cand_param_csv: string, file name.
        The column row in the example includes: solvent, D, P, H, epsilon, pi_star, HBD, HBA, logP_solv
        user can specify additional properties to be calculated while loading candidate database. 
        for example, hansen distance with starting material and the product
        kwargs:
        do_hansen_dist: bool
        hansen_dist_with: list of dict or one dict, keys are comp_name, D, P, H
        """
        all_solv_cand_df = pd.read_csv(all_cand_param_csv)
        #convert to dictionary
        self.all_cand_param_dict = all_solv_cand_df.to_dict(orient='records')
        self.all_cand_param_df = all_solv_cand_df
        print(f"Loading candidate parameters: success.\n")
        
    def _calc_hansen_dist(self, crt_df, hansen_dist_with_dict):
        """
        calculate hansen distance with current dataframe HSP columns
        hansen_dist_with_dict: dictionary
        """
        hansen_dist_with_name = hansen_dist_with_dict["comp_name"]
        hansen_dist_col_name = '_'.join(['r', hansen_dist_with_name])
        this_D = hansen_dist_with_dict["D"]
        this_P = hansen_dist_with_dict["P"]
        this_H = hansen_dist_with_dict["H"]
        
        this_r = np.sqrt(4*(crt_df["D"]-this_D)**2 + (crt_df["P"] - this_P)**2 + (crt_df["H"] - this_H)**2)
        
        crt_df[hansen_dist_col_name] = this_r
        
        print(f"Calculate hansen distance with {hansen_dist_with_name}: success.")
        return crt_df

    
    def load_exp_data(self, ip_exp_csv, **kwargs):
        """
        ip_exp_csv is the file name of input experimental data.
        The columns of the data file should follow the naming rule as rxn_id, solv_1, precent_1, solv_2, percent_2, Y
        There could be more co-solvents column; the name should follow solv_#, precent_# structure.
        Y is the observable column, in the example is the yield.
        Important rule: the name for solvent must match the name used in the all_cand_param data file.
        Mixtures parameters will be calculated and attached in each entry.
        """
        full_exp_df = pd.read_csv(ip_exp_csv)
        self.all_co_solv_col_name = self._get_co_solv_col_name(full_exp_df)
        full_exp_dict = full_exp_df.to_dict(orient="records")
        
        #operate each row/entry one by one
        #fetch parameters of corresponding solvent from database
        #then scale by percentage
        all_comb_exp = []
        
        for exp in full_exp_dict:
            
            rxn_id = exp["rxn_id"]
            print(f"Process {rxn_id} starts..")
            
            for i, solv_col in enumerate(self.all_co_solv_col_name["solv_col"]):
                this_solv_name = exp[solv_col]
                
                if type(this_solv_name) is str:
                
                    percent_col = self.all_co_solv_col_name["percent_col"][i]
                    this_solv_percent = exp[percent_col]

                    #fetch all parameters for this solvent
                    this_solv_params = self._fetch_solv_param(this_solv_name)["full_params"]
                    #apply this_solv_percent to all parameters
                    this_solv_scaled_params = self._scale_param_by_percent(this_solv_params, this_solv_percent)

                    if i == 0:
                        mix_solv_params = this_solv_scaled_params
                    else:
                        mix_solv_params = self._add_dict(this_solv_scaled_params, mix_solv_params)
            
            #process kwargs if any
            #calculate hansen distance
            if "do_hansen_dist" in kwargs:

                if kwargs["do_hansen_dist"]:

                    if type(kwargs["hansen_dist_with"]) is dict:
                        #calculate the hansen distance to the all_solv_cand_df
                        #this is the case when only one species need to be calculated
                        all_solv_cand_df = self._calc_hansen_dist(all_solv_cand_df, kwargs["hansen_dist_with"])
                        pass
                    elif type(kwargs["hansen_dist_with"]) is list:
                        #loop over the list and append calculated hansen distance to the all_solv_cand_df
                        #this is the case when multiple species need to be calculated
                        for with_this_comp in kwargs["hansen_dist_with"]:
                            all_solv_cand_df = self._calc_hansen_dist(mix_solv_params, with_this_comp)
            
            #include mix_solv_params to the experimental dataset
            comb_exp_dict = exp
            comb_exp_dict.update(mix_solv_params)
            
            all_comb_exp.append(comb_exp_dict)
        
        self.all_comb_exp = all_comb_exp
        print("\nLoading data: success.\n")
        
        #the meta data can be saved or not, user's choice
        if "save_meta_data" in kwargs:
            if kwargs["save_meta_data"]:
                meta_data_file_name = "op_exp_data_with_full_params.csv"
                meta_data_df = pd.DataFrame(self.all_comb_exp)
                self._save_meta_data(meta_data_df, meta_data_file_name, keep_index = False)
            
                
#         return self.all_comb_exp
    
    def _save_meta_data(self, meta_data_df, meta_data_file_name, keep_index = True):
        to_csv_path = self.op_folder_path + "\\" + meta_data_file_name
        meta_data_df.to_csv(to_csv_path, index = keep_index)
        print(f"{meta_data_file_name} saved.")
        
    
    
    def rm_zero(self):
        """
        remove zero yield data from column Y
        """
        non_zero_comb = []
        
        for entry in self.all_comb_exp:
            if entry['Y'] == 0:
                pass
            elif type(entry['Y']) not in [float, int]:
                err_rxn_id = entry["rxn_id"]
                raise ValueError(f"Wrong format found in the observable column. Please check entry: {err_rxn_id}. \n")
            else:
                non_zero_comb.append(entry)
        
        self.non_zero_comb = non_zero_comb
        print("Remove invalid/zero data: success.\n")
        self.non_zero_comb_df = pd.DataFrame(self.non_zero_comb)
        
        return self.non_zero_comb
                
        
    #parse column name and extract solv, precent part from the full columns
        
    def _get_co_solv_col_name(self, full_exp_df):
        """
        look at how many co-solvents are used in experiment.
        return column names for the solv/percent block
        """
        all_col = full_exp_df.columns
        #confirm how many co-solvents to parse
        all_co_solv = {
            "solv_col":[],
            "percent_col":[]
        }

        for i, col in enumerate(all_col):
            if "solv" in col.lower():
                this_solv_col = col
                this_percent_col = all_col[i+1]
                all_co_solv["solv_col"].append(this_solv_col)
                all_co_solv["percent_col"].append(this_percent_col)

        if len(all_co_solv["solv_col"]) == 0 or len(all_co_solv["percent_col"]) == 0:
            raise ValueError("Co-solvent or percent column was missing. Please check the input experimental data file.\n")
        elif len(all_co_solv["solv_col"]) != len(all_co_solv["percent_col"]):
            raise ValueError("Potential offset was detected in co-solvent and percentage info. Please check the input experimental data file.\n")
        
        #get total number of solvents in the mixture
        self.total_num_of_solv = len(all_co_solv["solv_col"])
        
        return all_co_solv

            
    
    def _scale_param_by_percent(self, full_param_dict, this_solv_percent):
        """
        scale each parameter by the percentage of current solvent
        """
        scaled_num_params = {}
        for item in full_param_dict.keys():
            this_value = full_param_dict[item]
            if type(this_value) is not str:
                scaled_num_params[item] = this_solv_percent * this_value
        return scaled_num_params
    
    def _add_dict(self, dict_1, dict_2):
        """
        dict_1 and dict_2 have same structure
        """
        mix_dict = {}
        for item in dict_1.keys():
            mix_dict[item] = dict_1[item] + dict_2[item]
        return mix_dict
        
    
    def _fetch_solv_param(self, solv_name):
        """
        fetch solvent parameters by solv_name
        """
        
        this_solv_params = {}
#         print(solv_name)
        
        for entry in self.all_cand_param_dict:
            #print(entry)
            
            if solv_name.lower() == entry["solvent"].lower():
                this_solv_params["this_solvent"] = solv_name
                this_solv_params["full_params"] = entry
        
#         print(this_solv_params)
        if not this_solv_params["this_solvent"]:
            raise ValueError("Target solvent absent from candidate database. Please check the solvent name or the database entry.")
        
        return this_solv_params
    
    
    #do corr plot, default Y with descriptors, can be modified by do any two columns
    def plot_corr(self, observable = "Y", plot_all = True, x_to_plot = ["D", "P", "H"], save_plot = True, do_std = False, **kwargs):
        """
        plot correlation between observable and a series of x
        if plot_all is true, default will plot every descriptor in the dataset,
        otherwise, only descriptors in x_to_plot will be plot.
        plot details can be specified in kwargs
        """
        #standardise the whole dataframe
        self.all_std_df = self._get_zsocre_for_full_df(self.non_zero_comb_df)# to be done in get_zscore function
        
        
        #if corr plot needs to be done with std data
        if do_std:
            full_data = self.all_std_df
        #or corr plot will be done by original data
        else:
            full_data = self.non_zero_comb_df
        
        if plot_all:
            #get all descriptors from full data list, find it after Y
            dscpt_col = self._get_full_descriptors_col()
            x_to_plot = dscpt_col
       
        for x_name in x_to_plot:
            sns.relplot(data = full_data, x = full_data[x_name], y = full_data[observable])
            plt.xlabel(x_name, fontsize = 15)
            plt.ylabel(observable, fontsize = 15)
            plt.tick_params(axis='both', which='major', labelsize=15)
            plt.tight_layout()
            
            #save all corr plots to the op_folder
            if save_plot:
                fig_to_save_name = self.op_folder_path + "\\corr_plot_" + x_name
                plt.savefig(fig_to_save_name)
        
        
        
        #do heatmap of descriptors
        dscpt_col = x_to_plot
        dscpt_col_df = full_data[dscpt_col]
        self.dscpt_df_corr = dscpt_col_df.corr()
        
        plt.clf()
        sns.heatmap(self.dscpt_df_corr, cmap = "coolwarm")
        
        if save_plot:
            heat_map_to_save_name = self.op_folder_path + "\\descriptor_heat_map"
            plt.tight_layout()
            plt.savefig(heat_map_to_save_name)

            print("All correlation plots have been saved.\n")
        
        #plot heatmap with observable included
        obs_with_dscpt = x_to_plot
#         print(observable)
#         print(obs_with_dscpt)
        obs_with_dscpt.append(observable)
#         print(obs_with_dscpt)

        obs_with_dscpt_df = full_data[obs_with_dscpt]
        obs_with_dscpt_df_corr = obs_with_dscpt_df.corr()
        
        plt.clf()
        sns.heatmap(obs_with_dscpt_df_corr, cmap = "coolwarm", square = True)
        
        
        if save_plot:
            plot_to_save_name = self.op_folder_path + "\\obs_with_descriptor_heat_map"
            plt.tight_layout()
            plt.savefig(plot_to_save_name)
            print("Correlation plot including observable has been saved. \n")
            plt.show()
            
        
        #save meta data
        if "save_meta_data" in kwargs:
            if kwargs["save_meta_data"]:
                #save corr matrix of descriptors
                self._save_meta_data(self.dscpt_df_corr, "corr_mat_of_descriptors_only.csv", keep_index = True)
                #save corr matrix of descriptors and observable
                self._save_meta_data(obs_with_dscpt_df_corr, "corr_mat_of_descriptors_with_observable.csv", keep_index = True)
        
        #pass the correlation matrix involving observable to the global scope
        self.obs_with_dscpt_df_corr = obs_with_dscpt_df_corr
        
        pass
    
    def _get_full_descriptors_col(self):
        """
        get column name of all descriptors, assumed to be after Y
        """
        all_col_name = self.non_zero_comb_df.columns
        dscpt_start_idx = 0
        for i, col_name in enumerate(all_col_name):
            if col_name == "Y":
                dscpt_start_idx = i+1
        
        dscpt_col = list(all_col_name[dscpt_start_idx:])
        #print(dscpt_col)
        if len(dscpt_col) != 0:
            print(f"Generate full descriptor: success.")
        
        self.num_of_dscpt = len(dscpt_col) #number of descriptors determine how many PC will be generated by default
        print(f"{self.num_of_dscpt} numerical descriptors are found. \n")
        
        return dscpt_col
    
    def _get_zsocre_for_full_df(self, ori_df, ignore_col = ["rxn_id", "percent_1", "percent_2", "percent_3", "percent_4"]):
        """
        calculate zscore for all the numerical columns in the original df
        """
        all_num_col = ori_df.select_dtypes(include = [np.number]).columns
        
        #process numerical columns but ignore percentage information
        to_zscore_col = []
        for col in all_num_col:
            if col not in ignore_col:
                to_zscore_col.append(col)
#         print(to_zscore_col)
                
        
        zscored_df = ori_df.copy()
        
        #store the zscore mean and std of experimental set
        #prepare for inverse pca
        zscore_info = {
            'col':[],
            'mean':[],
            'std':[]
        }
        
        
        #get descriptors zscored df
        dscpt_zscored_df = pd.DataFrame()
        dscpt_with_obs_zscored_df = pd.DataFrame()
        
        for c in zscored_df.columns:
            if c in to_zscore_col:
                
                #write zscore_info for statistical details
                zscore_info['col'].append(c)
                zscore_info['mean'].append(zscored_df[c].mean())
                zscore_info['std'].append(zscored_df[c].std())
                
                #write full zscored df
                zscored_df[c]= zscore(zscored_df[c])
                
                #write descriptors zsocred df with observable column
                dscpt_with_obs_zscored_df[c] = zscored_df[c]
                
                #write descriptors zscored df
                if c != "Y":
                    dscpt_zscored_df[c] = zscored_df[c]
        
        self.zscore_info = zscore_info
        
        self.dscpt_zscored_df = dscpt_zscored_df
        self.dscpt_with_obs_zscored_df = dscpt_with_obs_zscored_df
        
        return zscored_df
        
                
    def do_pca(self, n_comp = 0, to_plot = True, save_meta_data = True, **kwargs):
        """
        do pca on descriptors
        default not to do inverse pca
        if do inverse pca, pc_to_inv needs to be specified
        the dimension of pc to specified must match the n_comp
        n_comp default equals to the total number of descriptors
        when inverse pca needs to be carried out,
        n_comp needs to be specified as selected number of pc
        """
        #confirm number of PCs
        if not n_comp:
            n_comp = self.num_of_dscpt
        #print(n_comp)
        
        
        self.pca_log = {}
        self.pca_log['n_comp'] = n_comp

        #prepare PCA input data
        #PCA input data consist of zscored descriptors data
        pca_ip_data = self.dscpt_zscored_df
        self.pca_log['pca_ip_data'] = pca_ip_data
        self.pca_log['descriptors'] = list(pca_ip_data.columns)

        #status info
        print(f"PCA starts...")
        print(f"Original descriptors:{list(pca_ip_data.columns)}")

        #do pca
        pca = PCA(n_components = n_comp)
        pca.fit(pca_ip_data)
        self.pca_log['pca'] = pca

        pca_comp = pca.components_
        self.pca_log['pca_comp'] = pca_comp

        #explained variance ratio
        exp_var_ratio = pca.explained_variance_ratio_
        self.pca_log['exp_var_ratio'] = exp_var_ratio
        #print(f"Explained variance ratio: {exp_var_ratio}")

        #transform
        transformed_data = pca.fit_transform(pca_ip_data)
        self.pca_log['trans_data'] = transformed_data
        #print(transformed_data)

        #loading_matrix
        loadings = pca.components_.T * np.sqrt(pca.explained_variance_)
        self.pca_log['loadings'] = loadings

        #prepare column name after PC
        loading_matrix_column = []
        for i in range(0, n_comp):
            loading_matrix_column.append('PC'+str(i+1))

        self.pca_log['loading_mat_col'] = loading_matrix_column

        #create loading matrix
        loading_matrix = pd.DataFrame(loadings, columns = loading_matrix_column, index = pca_ip_data.columns)
        #print(loading_matrix)
        self.pca_log['loading_mat'] = loading_matrix

        #get explained_variance
        exp_var = pca.explained_variance_
        #print(f"explained variance: {exp_var}")
        self.pca_log['exp_variance'] = exp_var

        #get score or transformed data
        X_pca = pca.transform(pca_ip_data)
        scores = X_pca
        #print(f"score: {scores}")
        self.pca_log['scores'] = scores


        #print status info
        print(f"PCA completed. \n Please refer to pca_log for details.")

        #plot things and save or not...
        if to_plot:
            plt.clf()

            #scree plot of explained variance ratio versus PC
            plot_name = "pca_scree_plot_n_" + str(n_comp)
            print(f"{plot_name}:\n")
#             to_plot_dict = dict(map(lambda i, j : (i, j), self.pca_log['loading_mat_col'], self.pca_log['exp_var_ratio']))

            to_plot_data = pd.DataFrame()
            to_plot_data['PC'] = self.pca_log['loading_mat_col']
            to_plot_data['exp_var_ratio'] = self.pca_log['exp_var_ratio']
#             print(to_plot_data)

            fig, ax = plt.subplots()
            ax = sns.barplot(x = to_plot_data['PC'], y = to_plot_data['exp_var_ratio'])
            plt.xlabel('PC')
            plt.ylabel('Explained variance ratio')
            plot_to_save_as = self.op_folder_path + "\\"+ plot_name
            plt.savefig(plot_to_save_as)
            plt.show()

            if save_meta_data:
                file_name = plot_name + ".csv"
                self._save_meta_data(to_plot_data, file_name, keep_index = False)


            #plot loading matrix
            plt.clf()
            plot_name = "pca_loading_matrix_n_" + str(n_comp)
            print(f"{plot_name}:\n")

            fig, ax = plt.subplots(figsize = (10, 10))
            ax = sns.heatmap(self.pca_log['loading_mat'], annot = True, cmap = 'coolwarm')
            plot_to_save_as = self.op_folder_path + "\\" + plot_name
            plt.savefig(plot_to_save_as)
            plt.show()

            if save_meta_data:
                file_name = plot_name + ".csv"
                self._save_meta_data(self.pca_log['loading_mat'], file_name, keep_index = True)

            
            

    def get_pc_map(self, n_comp = 2, with_y = True, save_plot = True, save_meta_data = True):
        """
        plot observable in selected PC space, default first two PC
        """
        all_scores = self.pca_log['scores'] #nparray
        all_loadings = np.array(self.pca_log['loading_mat']) #dataframe to array
        
        if n_comp < 2:
            raise ValueError("Number of PC must be equal or over 2.")
        elif n_comp >= len(all_scores):
            raise ValueError(f"Number of PC must be less than {len(all_scores)}.")
        else:
            total_pc = n_comp
            plot_df = pd.DataFrame(all_scores, columns = self.pca_log['loading_mat_col'])
            plot_df['Y'] = self.all_std_df["Y"]
            
            #plot consecutive two PC space
            #max column index of PC is n_comp minus 1
            #for example, if n_comp equals 3, then two 2D spaces: col0-col1, col1-col2 needs to be created
            for pc_col in range(n_comp-1):
                
                fig, ax = plt.subplots(figsize = (8, 8))
                
                x_pc_id = pc_col #pc_x
                y_pc_id = pc_col + 1 #pc_y
                
                #get score
                x_score = all_scores[:, x_pc_id]
                y_score = all_scores[:, y_pc_id]
                
                #plot score
                if not with_y:
                    plt.scatter(x_score, y_score)
                else:
                    col_x = "PC" + str(x_pc_id + 1)
                    col_y = "PC" + str(y_pc_id + 1)
#                     print(plot_df.columns)
#                     print(plot_df[col_y])
                    sns.scatterplot(data = plot_df, x = col_x, y = col_y, hue = plot_df["Y"], hue_norm = (0,1), s = 100, palette = "crest")
#                     ax = sns.scatterplot(plot_df, x = col_x, y = 'PC2', hue = 'Y', s = 100, palette = "crest")
                    ax.legend(bbox_to_anchor = (1.2, 1))
                    plt.title('Standardised yield')
                
                #plot the loadings
                for i, v in enumerate(all_loadings):
                    ax.arrow(0, 0, v[x_pc_id], v[y_pc_id], head_width=0.05, head_length=0.05)
                    ax.text(v[x_pc_id] * 1.2, v[y_pc_id] * 1.2, self.pca_log['descriptors'][i])
                
                ax.set_xlabel(f"PC{x_pc_id + 1}", fontsize = 20)
                ax.set_ylabel(f"PC{y_pc_id + 1}", fontsize = 20)
                plt.xticks(fontsize = 15)
                plt.yticks(fontsize = 15)
                
                #save biplot
                if save_plot:
                    plt.tight_layout()
                    if not with_y:
                        plot_name = f"pca_map_{n_comp}d_biplot_PC{x_pc_id + 1}_PC{y_pc_id + 1}"
                    else:
                        plot_name = f"pca_map_with_obs_{n_comp}d_biplot_PC{x_pc_id + 1}_PC{y_pc_id + 1}"
                    plot_to_save_as = self.op_folder_path + "\\" + plot_name
                    plt.savefig(plot_to_save_as)
                    plt.show()
                
                #plot zoomed-in loadings for descriptors
                #i.e., loading only, but large
                
                fig, ax = plt.subplots(figsize = (8, 8))
                for i, v in enumerate(all_loadings):
                    ax.arrow(0, 0, v[x_pc_id], v[y_pc_id], head_width=0.02, head_length=0.05, facecolor = "black")
                    ax.text(v[x_pc_id] * 1.15, v[y_pc_id] * 1.15, self.pca_log['descriptors'][i], fontsize = 20)
                ax.set_xlabel(f"PC{x_pc_id + 1}", fontsize = 20)
                ax.set_ylabel(f"PC{y_pc_id + 1}", fontsize = 20)
                ax.set_xlim(-1.3, 1.3)
                ax.set_ylim(-1.3, 1.3)
                plt.xticks(fontsize = 15)
                plt.yticks(fontsize = 15)
                plt.tight_layout()
                
                #save loading plot
                if save_plot:
                    plot_name = f"pca_map_{n_comp}d_loading_PC{x_pc_id + 1}_PC{y_pc_id + 1}"
                    plot_to_save_as = self.op_folder_path + "\\" + plot_name
                    plt.savefig(plot_to_save_as)
                    plt.show()
                
    def inv_pca_coord(self, n_comp = 2, pc_to_inv = [5, -1], save_meta_data = True, **kwargs):

        #inverse pca

        print("Inverse PCA starts...")
        if len(pc_to_inv) == n_comp:
            #calculate inverse pca for all descriptors

            #convert input coordinate to np.array
            pc_to_inv = np.array(pc_to_inv)
            
            
           #generate transform matrix for inverse pca set
            pca_for_inv = PCA(n_components = n_comp)
            pca_for_inv.fit(self.dscpt_zscored_df)

            #calculate std descriptors in the original high dim space
            inv_std_dscpts = pca_for_inv.inverse_transform(pc_to_inv)

            print(f"Inverse PCA done.")
            print(f"{self.pca_log['descriptors']}\n for PC coordinate: {pc_to_inv} in the original dimension space is:\n {inv_std_dscpts}.\n")
            
            inv_pca_log = {}
            
#             self.inv_pca_log['pc_to_inv'] = pc_to_inv
#                 self.pca_log['inv_pca_col']
            inv_pca_log['descriptors'] = self.pca_log['descriptors']
            inv_pca_log['inv_pca_orig_space_value'] = inv_std_dscpts
            
#             print(inv_pca_log)
            inv_pca_log_df = pd.DataFrame.from_dict(inv_pca_log)
            
            
            if save_meta_data:
                
                pc_to_inv_list = list(pc_to_inv)
                coord_in_name = '_'.join(str(coord) for coord in pc_to_inv_list)
                    
                file_name = f"inverse_pca_df_{coord_in_name}.csv"
                self._save_meta_data(inv_pca_log_df, file_name, keep_index = False)
            
            self.inv_pca_full_log = {}
            self.inv_pca_full_log['n_comp'] = n_comp
            self.inv_pca_full_log['coord_to_inv'] = pc_to_inv
            self.inv_pca_full_log['descriptors'] = inv_pca_log['descriptors']
            self.inv_pca_full_log['orig_space_value'] = inv_std_dscpts
                


        else:
            raise ValueError("The dimension of PC coordinate to inverse must match the number of PC.")

                    
    def sel_feature(self, n_comp = 2, pc_to_rank = 'PC1', abs_value = True):
        """
        rank the feature by PCA loadings
        total PC number should match the n_comp in inv_pca_coord
        if abs_value is True, feature will be ranked by absolute value of projection
        """
        #do pca for rank purpose
        pca_for_rank = PCA(n_components = n_comp)
        pca_for_rank.fit(self.dscpt_zscored_df)
        
        #loading on pc_to_rank is required
        #first get loading matrix
        loadings = pca_for_rank.components_.T * np.sqrt(pca_for_rank.explained_variance_)
        #prepare column name after PC
        loading_matrix_for_rank_column = []
        for i in range(0, n_comp):
            loading_matrix_for_rank_column.append('PC'+str(i+1))
        #create loading matrix
        loading_matrix_for_rank = pd.DataFrame(loadings, columns = loading_matrix_for_rank_column, index = self.dscpt_zscored_df.columns)
        
        #select the PC to rank
        if pc_to_rank in loading_matrix_for_rank.columns:
            loading_col_to_rank = loading_matrix_for_rank[pc_to_rank]
            #sort from highest to lowest
            if not abs_value:
                #feature will be ranked by absolute value of projection
                desc_loading_col = loading_col_to_rank.sort_values(ascending = False)
                print(f"The projection of all descriptors on {pc_to_rank} axis is: \n")
            
            else:
                abs_loading_col_to_rank = abs(loading_col_to_rank)
                desc_loading_col = abs_loading_col_to_rank.sort_values(ascending = False)
                print(f"The absolute projection of all descriptors on {pc_to_rank} axis is: \n")
            
            
            print(desc_loading_col)

            self.sel_feature_log = {}
            self.sel_feature_log["pc_to_rank"] = pc_to_rank
            self.sel_feature_log["do_abs"] = abs_value
            self.sel_feature_log["loading_col_to_rank"] = loading_col_to_rank
            self.sel_feature_log["ranked_loading_col"] = desc_loading_col
            
            
            #point out the highest projection
            best_dscpt = desc_loading_col.index[0]
            best_dscpt_value = desc_loading_col[0]
            
            print(f"The best descriptor is: {best_dscpt}. \n Projection: {best_dscpt_value}. \n")
            
            self.sel_feature_log["best_descriptor"] = best_dscpt
            self.sel_feature_log["best_descriptor_projection"] =  best_dscpt_value
            
            print("The full feature selection log can be found by calling: sel_feature_log\n")
        else:
            raise ValueError("pc_to_rank not found in loading matrix.")
        
    def hop_back_to_high_dim(self, sel_by = 0):
        """
        By the time this function is called, we assume inverse pca and feature selection has both been completed.
        sel_by will be passed with the best descriptor determined in the feature selection step,
        whereas user can specify feature name, such as "r_sm", if in need
        The output will tell you what is the standardised value of the selected feature in the original high dimensional space
        """
        if self.sel_feature_log and self.inv_pca_full_log:
            if sel_by == 0:
                #default taking the best descriptor recorded in the feature selection step
                sel_by = self.sel_feature_log["best_descriptor"]
            else:
                if sel_by not in self.dscpt_zscored_df.columns:
                    raise ValueError("Please choose a valid descriptor.")
            
            #convert full inverse pc log to dataframe:
            inv_pca_dscpt_ori_space_dict = {}
            inv_pca_dscpt_ori_space_dict["descriptors"] = self.inv_pca_full_log['descriptors']
            inv_pca_dscpt_ori_space_dict["orig_space_value"] = self.inv_pca_full_log['orig_space_value']
            
            inv_pca_dscpt_ori_space_df = pd.DataFrame.from_dict(inv_pca_dscpt_ori_space_dict)
            inv_pca_dscpt_ori_space_df = inv_pca_dscpt_ori_space_df.set_index("descriptors")
            
            print(inv_pca_dscpt_ori_space_df)
            
            orig_feature_value = inv_pca_dscpt_ori_space_df[inv_pca_dscpt_ori_space_df.index == sel_by]["orig_space_value"][0]
            
#             print(orig_feature_value)
            
            print(f"The standardised value of {sel_by} in the original feature space is: {orig_feature_value}.\n")
            
            #write to log
            self.high_dim_log = {}
            self.high_dim_log["sel_by"] = sel_by
            self.high_dim_log["orig_feature_value"] = orig_feature_value
            
            print(f"Results have been stored in high_dim_log. \n")
            
            
        else:
            raise ValueError("Please make sure inverse pca and feature selection has been done")
    
    #do solvent prediction based on target selected feature
    
    #by default, features that deviates least from the target parameter will be selected as candidates
    #alternative method is to specify two solvents and create the mixture
    #a more inclusive option is to use the full candidate list, i.e., solv_pred, but not hsp
    def pred_mixture(self, n_solv = 2, target_feature = "default", target_value = "default", pred_mode = "closest", save_meta_data = True, **kwargs):
        """
        n_solv: number of solvents to be included in the mixture, default 2
        target_feature: default self.high_dim_log["sel_by"]
        target_value: default self.high_dim_log["orig_feature_value"]
        pred_mode: default "closest" -> find the closest candidates based on the difference between std_cand_set from the target
        pred_mode: "specify" -> specify target candidates using keywords: solv_cand = [solv_cand_name_1, solv_cand_name_2, ...]
        pred_mode: "full" -> prediction will be made by all possibilities of candidate combinations.
        save_meta_data: default True, predicted solvent combinations and calculated feature value and error information will be stored in two spreadsheets.
        """
        
        if target_feature == "default":
            target_feature = self.high_dim_log["sel_by"]
        elif target_feature not in self.dscpt_zscored_df.columns:
            raise ValueError("Invalid target feature.")
        
        if target_value == "default":
            target_value = self.high_dim_log["orig_feature_value"]
        elif type(target_value) not in [int, float]:
            raise ValueError("Invalid target value.")
        
        #prepare std_cand_set
        self._scale_cand_set()
        
        
        if pred_mode == "closest":
            n_cands, n_values = self._find_closest_n(n_solv, target_feature, target_value)
            print(n_cands, n_values)
            
            self._calc_norm_solv_coeff(n_solv, n_values, target_value)

            #map cand with ratio
            solv_coeff_dict = self._gen_solv_coeff_dict(n_cands)
            

            if save_meta_data:
                solv_coeff_df = pd.DataFrame.from_dict(solv_coeff_dict)
                solv_coeff_file_name = f"mixture_component_for_{target_feature}_{target_value:.2f}_{n_solv}_cand.csv"
                solv_coeff_calc_detail_df = pd.DataFrame.from_dict(self.solv_coeff_calc_detail)
                solv_calc_detail_file_name = f"mixture_calculation_info_for_{target_feature}_{target_value:.2f}_{n_solv}_cand.csv"
                self._save_meta_data(solv_coeff_df, solv_coeff_file_name, keep_index = False)
                self._save_meta_data(solv_coeff_calc_detail_df, solv_calc_detail_file_name, keep_index = False)

                
                

            
        
        elif pred_mode == "specify":
            #specify a list of solvent candidate in kwargs
            if "solv_cand" in kwargs:
                solv_cand = kwargs["solv_cand"]
            
                if len(solv_cand) != n_solv:
                    raise ValueError(f"Please ensure the number of proposed candidates equal to the number of solvents: {n_solv}. \n")
                
                #extract feature value for all solv_cand provided
                n_cands, n_values = self._get_value_for_solvent(solv_cand, target_feature)
                print(n_cands, n_values)

                self._calc_norm_solv_coeff(n_solv, n_values, target_value)
                #get error

                #map cand with ratio
                solv_coeff_dict = self._gen_solv_coeff_dict(n_cands)


            if save_meta_data:
                solv_coeff_df = pd.DataFrame.from_dict(solv_coeff_dict)
                solv_coeff_file_name = f"mixture_component_for_{target_feature}_{target_value:.2f}_{n_solv}_cand.csv"
                solv_coeff_calc_detail_df = pd.DataFrame.from_dict(self.solv_coeff_calc_detail)
                solv_calc_detail_file_name = f"mixture_calculation_info_for_{target_feature}_{target_value:.2f}_{n_solv}_cand.csv"
                self._save_meta_data(solv_coeff_df, solv_coeff_file_name, keep_index = False)
                self._save_meta_data(solv_coeff_calc_detail_df, solv_calc_detail_file_name, keep_index = False)

            
            else:
                raise ValueError(f"Please specify {n_solv} solvent candidates using keyword: solv_cand = [your solvent names]. \n")
        
        elif pred_mode == "full":
            #do full prediction by iterating through all candidates
            #see if we can connect to solv_pred
            self._calc_full_cand_mixture(n_solv, target_feature, target_value)
            if self.full_mixture_dict:
                
                if save_meta_data:
                    full_mixture_df = pd.DataFrame.from_dict(self.full_mixture_dict)
                    full_mixtue_file_name = f"mixture_component_for_{target_feature}_{target_value:.2f}_{n_solv}_all_cand.csv"
                    self._save_meta_data(full_mixture_df, full_mixtue_file_name, keep_index = False)

                    #get valid results
                    vld_full_mixture_df = full_mixture_df.loc[full_mixture_df['valid'] == "True"]
                    vld_full_mixture_file_name = f"mixture_component_for_{target_feature}_{target_value:.2f}_{n_solv}_all_cand_valid.csv"
                    self._save_meta_data(vld_full_mixture_df, vld_full_mixture_file_name, keep_index = False)

                    if len(vld_full_mixture_df['group']) == 0:
                        print("Warning: No valid combinations can be predicted.\n Please try a different target. \n") 
            
        
        else:
            raise ValueError("Invalid pred_mode. Options: [closest, specify, full]. \n")

        pass
    
    #copy the candidate dataframe and scale by experimental zscore parameters
    def _scale_cand_set(self):
        """
        make a copy of candidate set and convert all values by same parameters used in experimental data zscored
        The aim is to make candidate solvent parameters comparable to the information processed by PCA
        """
        #extract the zscore info
        all_feature_name = self.zscore_info["col"]
        all_mean = self.zscore_info["mean"]
        all_std = self.zscore_info["std"]
        
        #revisit the candidate set
        std_cand_df = self.all_cand_param_df.copy()
        
        for i, feature in enumerate(all_feature_name):
            
            if feature in std_cand_df.columns:
                this_mean = all_mean[i]
                this_std = all_std[i]
                orig_feature_value = std_cand_df[feature]
                zscore_feature_value = (orig_feature_value - this_mean)/this_std
                std_cand_df[feature] = zscore_feature_value

        self.std_cand_set = std_cand_df.set_index("solvent")
        
        
        print("Standardisation on candidate set based on zscore parameters used for experimental set: done. \n")
        print("Standardised candidate set: \n")
        print(self.std_cand_set)
        print("The standardised candidate set can be found in std_cand_set. \n")
        
    
    def _find_closest_n(self, n_solv, target_feature, target_value):
        """
        find the closest n solvents from the std_cand_set from the target feature value
        """
        #calculate the absolute difference between target feature value and std_cand_set
        feature_col_in_cand_set = self.std_cand_set[target_feature]
        diff_of_cand_feature_from_target = abs(feature_col_in_cand_set - target_value)
        #sort by ascending
        asc_diff = diff_of_cand_feature_from_target.sort_values(ascending = True)
        
        print(f"The absolute difference of {target_feature} in the original dataset from {target_value} is: \n")
        print(asc_diff)
        
        first_n_row = asc_diff.iloc[:n_solv]
        print(f"\n First {n_solv} closest candidates: \n")
        print(first_n_row)
        first_n_cands = list(first_n_row.index)
        n_cands, first_n_values = self._get_value_for_solvent(first_n_cands, target_feature)
        
#         print(first_n_cands, first_n_value)
        
        return first_n_cands, first_n_values
    
    def _get_value_for_solvent(self, solv_name, target_feature):
        """
        find feature value for solvent
        solv_name is a list
        """
#         feature_col_in_cand_set = self.std_cand_set[target_feature]
        n_cands = []
        n_values = []
        
        for solv in solv_name:
            if solv in self.std_cand_set.index:
                n_cands.append(solv)
                n_values.append(self.std_cand_set[self.std_cand_set.index == solv][target_feature][0])
            else:
                raise ValueError(f"{solv} is not in candidate set. \n")
        
        return n_cands, n_values
    
    
    def _calc_solv_coeff(self, n_solv, n_values, target_value):
        """
        work out coefficients of each solvent
        """
        feature_mat = np.array(n_values)
        reshaped_feature_mat = feature_mat.reshape(1, n_solv)
        feature_mat_pinv = np.linalg.pinv(reshaped_feature_mat)
        coeff_solv = np.dot(feature_mat_pinv, target_value)

        # print(coeff_solv)
        if all(coeff_solv) >= 0:
            norm_coeff_solv = self._norm_coeff(coeff_solv)
            return norm_coeff_solv
        
        else:
            #invalid result, something is less than 0
            return
        

    def _norm_coeff(self, all_coeff):
        """
        normalise each coefficient
        """     
        total_coeff = sum(all_coeff)
        norm_coeff = []
        for coeff in all_coeff:
            new_coeff = coeff/total_coeff
            norm_coeff.append(new_coeff)
        return norm_coeff

    def _calc_norm_solv_coeff(self, n_solv, n_values, target_value):
            
        print("Calculating mixture coefficient...\n")
        self.all_solv_coeff = self._calc_solv_coeff(n_solv, n_values, target_value)


        if self.all_solv_coeff:
            print("Done.\n")
            print(f"Coefficient matrix: {self.all_solv_coeff}.\n")

            real_rst, calc_error, calc_error_percent = self._calc_error(n_values, target_value)
            self.solv_coeff_calc_detail = {
                'real_rsl': real_rst,
                'calc_error' : calc_error,
                'percent_error' : calc_error_percent
            }
            print("Calculation details: \n")
            print(self.solv_coeff_calc_detail)


        else:
            print("No valid prediction can be made.\n")

    def _calc_error(self, n_values, target_value):
        """
        calculate error of predicted results with target value
        """
        real_rst = np.dot(np.array(n_values), self.all_solv_coeff)
        calc_error = real_rst - target_value
        calc_error_percent = calc_error/target_value
        return real_rst, calc_error, calc_error_percent
        

    def _gen_solv_coeff_dict(self, n_cands):
        """
        map solvent candidate name with as-calculated ratio
        """

        solv_coeff_dict = {
            'cand' : [],
            'ratio' : []
            }
        
        for i, cand in enumerate(n_cands):
            solv_coeff_dict['cand'].append(cand)
            solv_coeff_dict['ratio'].append(self.all_solv_coeff[i][0])
        
        print(solv_coeff_dict)
        return solv_coeff_dict
    

    def _get_all_solv_cand_comb(self, n_solv):
        """
        iterate through the candidate set and pick out all n_solv combinations
        """
        all_solv_cand = list(self.std_cand_set.index) # a list of all candidate solvents
        # print(all_solv_cand)
        all_solv_comb = combinations(all_solv_cand, n_solv)
        return all_solv_comb
    
    def _calc_full_cand_mixture(self, n_solv, target_feature, target_value):
        """
        calculate solvent mixtures based on all combinations for the "full" pred_mode
        """
        all_solv_comb = self._get_all_solv_cand_comb(n_solv)
        
        #iterate and do normal mixture calculation cycle
        
        #initialise result dict
        full_mixture_dict = {
            "group":[]
        }

        for n in range(n_solv):
            solv_col_name = f"solvent_{n+1}"
            full_mixture_dict[solv_col_name] = []
            percent_col_name = f"percent_{n+1}"
            full_mixture_dict[percent_col_name] = []
        
        full_mixture_dict["real_result"] = []
        full_mixture_dict["error"] = []
        full_mixture_dict["error_percent"] = []
        full_mixture_dict["valid"] = []

        for i, solv_comb in enumerate(all_solv_comb):

            solv_name = list(solv_comb) #convert tuple to list
            n_cands, n_values = self._get_value_for_solvent(solv_name, target_feature) #get cand name and original feature value
            self._calc_norm_solv_coeff(n_solv, n_values, target_value)
            solv_coeff_dict = self._gen_solv_coeff_dict(n_cands) # coeff of this combination detail
            solv_coeff_calc_detail_dict = self.solv_coeff_calc_detail #error and other calc detail

            if solv_coeff_dict and solv_coeff_calc_detail_dict:
                # full_mixture_dict["valid"].append("True")

                perc_vld = True


                for n in range(n_solv):
                    solv_n_name = solv_coeff_dict['cand'][n]
                    solv_n_ratio = solv_coeff_dict['ratio'][n]

                    if solv_n_ratio < 0:
                        perc_vld = False

                    solv_col_name = f"solvent_{n+1}"
                    full_mixture_dict[solv_col_name].append(solv_n_name)
                    percent_col_name = f"percent_{n+1}"
                    full_mixture_dict[percent_col_name].append(solv_n_ratio)

                    
                
                real_result = solv_coeff_calc_detail_dict["real_rsl"][0]
                error = solv_coeff_calc_detail_dict["calc_error"][0]
                error_percent = solv_coeff_calc_detail_dict["percent_error"][0]

                full_mixture_dict["group"].append(i+1) #write group number
                full_mixture_dict["real_result"].append(real_result)
                full_mixture_dict["error"].append(error)
                full_mixture_dict["error_percent"].append(error_percent)
                full_mixture_dict["valid"].append(str(perc_vld))
            
            
            else:
                pass
                # full_mixture_dict["valid"].append("False")
        if len(full_mixture_dict["group"]) != 0:
            print("Full prediction: done.\n")
            self.full_mixture_dict = full_mixture_dict
            # return full_mixture_dict
        else:
            raise ValueError("No valid combinations. Please try another target feature value. \n")
