"""
{Test script for scPROFILE v04 - with NBU scRNA-seq dataset}
"""
__author__ = 'Saran PANKAEW'
__version__ = '0.4.5'
__maintainer__ = 'Saran PANKAEW'
__email__ = 'saran.pankeaw@curie.fr'
__status__ = 'development'
__date__ = '04/04/2023'

###########
#TO DO ####

# Check how to process mutation profile
# Check if the data processing from R and python gives the same result

###########


import pandas as pd
import maboss
from sklearn import mixture
from scipy import stats
from  scipy.stats import kurtosis, median_abs_deviation
# from unidip.dip import diptst 
import diptest as dptest
import numpy as np
from tqdm import tqdm
import multiprocessing as mp
import time
import os
#import shutil
from datetime import datetime


def percent_0s(df, axis=0):
    n_0s = np.sum(df==0, axis=axis)
    return n_0s/df.shape[axis]

def discard_percent_0s(df, max_percent=0.95,axis=0):
    p_0s = percent_0s(df, axis=axis)
    if axis==0:
        df_filt = df.transpose()[p_0s<=max_percent]
        return df_filt.transpose()
    else:
        return df[p_0s<=max_percent]
    
#function combinating the three methods to define binarizable variables used in PROFILE
def binarizable(df, bi_min=None, kurtosis_max=None, diptest_max=None, return_df=False):
    df_ = df.copy()
    if (bi_min==None) & (kurtosis_max==None) & (diptest_max==None):
        print('No filter will be applied since no threshold has been precised.')
        return None
    if bi_min != None:
        df_ = df_.transpose()[df_.bimodality_index()>=bi_min].transpose()
    if kurtosis_max != None:
        df_ = df_.transpose()[kurtosis(df_)<kurtosis_max].transpose()
    if diptest_max!=None:
        df_ = df_.transpose()[diptest(df_)<diptest_max].transpose()
    if return_df==False:
        return df_.columns
    else:
        return df_

#method to calculate the bimodality index for all the variables of a dataframe
def bimodality_index(df, axis=0):
    #don't forget that mixture gives a dofferent result at each run since the process starts from a random state.
    BIs = []
    if axis==0:
        for i in range(df.shape[1]):            
            g = mixture.GaussianMixture(n_components=2)
            g.fit(np.array(df.iloc[:,i]).reshape(-1, 1))
            
            sigma = np.sqrt(np.mean(g.covariances_))
            delta = abs(g.means_[0] - g.means_[1])/sigma
            pi = g.weights_[0]
            BI = delta * np.sqrt(pi*(1-pi))
            BIs.append(BI)
    else:
        for i in range(df.shape[0]):        
            g = mixture.GaussianMixture(n_components=2)
            g.fit(np.array(df.iloc[i,:]).reshape(-1, 1))
                  
            sigma = np.sqrt(np.mean(g.covariances_))
            delta = abs(g.means_[0] - g.means_[1])/sigma
            pi = g.weights_[0]
            BI = delta * np.sqrt(pi*(1-pi))
            BIs.append(BI)
    return np.array(BIs)
#add the function modality_index as a method of pandas.DataFrame objects.
pd.DataFrame.bimodality_index = bimodality_index

def diptest(df, axis=0):
    dip_pvalue = []
    if axis==0:
        print(df.shape)
        for i in range(df.shape[1]):
            gene = df.iloc[:,i]
        
            # dip_pvalue.append(diptst(gene)[1])
            dip_pvalue.append(dptest.diptest(gene)[1])
        dip_pvalue = pd.Series(dip_pvalue, index=df.columns)
        return dip_pvalue
    else:
        raise
            

def binarize(df, axis=0):
    bin_df = df.transpose()[[False for n in range(df.shape[1])]].transpose()
    if axis==0:
        for i in range(df.shape[1]):
            gene = df.iloc[:,i]
            g = mixture.GaussianMixture(n_components=2)
            g.fit(np.array(df.iloc[:,i]).reshape(-1, 1))
            if g.means_[0]<g.means_[1]:
                #print(1)
                bin_df[df.columns[i]]=g.predict(np.array(df.iloc[:,i]).reshape(-1, 1))
            else:
                #print(2)
                bin_df[df.columns[i]] = np.abs(g.predict(np.array(df.iloc[:,i]).reshape(-1, 1))-1)

    else:
        bin_df = bin_df.transpose()
        for i in range(df.shape[0]):
            gene = df.iloc[i,:]
            
            g = mixture.GaussianMixture(n_components=2)
            g.fit(np.array(df.iloc[:,1]).reshape(-1, 1))
            g.predict(np.array(df.iloc[:,1]).reshape(-1, 1))
            bin_df[df.index[i]]=g.predict(np.array(df.iloc[:,i]).reshape(-1, 1))
    return bin_df


def uni_norm(df):
    df_ = df.copy()

    lambd = np.log(3)/median_abs_deviation(df_)
    df_norm = 1/(1+np.exp(-lambd*(df_-df_.median())))
    return df_norm
def inflated0_norm(df):
    if len(df.columns)!=0:
        percent_1 = np.percentile(df, 1)
        percent_99 = np.percentile(df, 99)
        df_norm = (df - percent_1)/(percent_99 - percent_1)
        df_norm[df_norm>1]=1
        df_norm[df_norm<0]=0
    else:
        df_norm=df
    return df_norm
def inflated0_test(df):
    amplitudes = np.max(df) - np.min(df)
    peaks=[]
    results = []    
    for i in range(df.shape[1]):
        gene = df.iloc[:,i]
        kde1 = stats.gaussian_kde(gene)
        
        x_eval = np.linspace(np.min(gene), np.max(gene), 1000)
        proba = kde1(x_eval)
        
        peaks.append(x_eval[np.where(proba==max(proba))[0]])
        if peaks[-1]<amplitudes[i]/10:
            results.append(True)
        else:
            results.append(False)
    return df.columns[results]

######################################################################
##                 Personnalized data calculation                   ##
######################################################################

def mutations_effects_compilation(fnames_celllines_mutations, fname_oncoKB_database, homozygous_loss_only=True):
    oncoKB = pd.read_csv(fname_oncoKB_database, sep='\t', index_col=0)
    
    
    
    cosmic_loss = ['Substitution - Nonsense', 'Deletion - Frameshift',\
                   'Insertion - Frameshift', 'Complex - frameshift']
    cosmic_unknown = ['Unknown', 'Substitution - coding silent', 'Nonstop extension', 'Deletion - In frame', 'Insertion - In frame']
    #keep In frame insertion/deletion in unknown or check somewhere? 
    oncoKB_loss = ['Likely Loss-of-function', 'Loss-of-function']
    oncoKB_gain = ['Likely Gain-of-function', 'Gain-of-function']
    oncoKB_switch = ['Switch-of-function', 'Likely Switch-of-function']
    oncoKB_neutral = ['Neutral']
    oncoKB_undefined = [oncoKB['Mutation Effect'].unique()[9],'Inconclusive']

    cellline_names = ['Car1', 'HT29', 'LS411N', 'SW1417', 'SW1463', 'SW403', 'SW480', 'SW620', 'SW837', 'HCT116']
    for cellline_name in fnames_celllines_mutations:

        mutations = pd.read_csv(fnames_celllines_mutations[cellline_name], index_col = 0)
        mutations.index = [ind.split('_')[0] for ind in mutations.index]



        functional_impact = []
        for i in range(len(mutations['AA Mutation'])):
            mutation = mutations.iloc[i,:]
            gene_name = mutation.name
            mutation_name = mutation['AA Mutation'].split('.')[1]

            if mutation['Type'] in cosmic_loss:
                predicCOSMIC = 'OFF'
                functional_impact.append(predicCOSMIC)

            elif mutation['Type']=='Substitution - Missense': 
                try:
                    oncoKB_mutation = oncoKB.loc[gene_name][oncoKB.loc[gene_name]['Alteration'] == mutation_name]
                    oncoKB_impact = oncoKB_mutation['Mutation Effect'][0]
                    if oncoKB_impact in oncoKB_gain:
                        predicKB = 'ON'
                    elif oncoKB_impact in oncoKB_loss:
                        predicKB = 'OFF'
                    elif oncoKB_impact in oncoKB_neutral:
                        predicKB = 1/2
                    elif oncoKB_impact in oncoKB_switch:
                        predicKB = 2
                    elif oncoKB_impact in oncoKB_undefined:
                        predicKB = -1
                except:
                    predicKB = -12
                functional_impact.append(predicKB)

            elif mutation['Type'] in cosmic_unknown:
                predicCOSMIC = -13
                functional_impact.append(predicCOSMIC)


        mutations['Impact'] = functional_impact
        gain = mutations[mutations['Impact']=='ON']
        loss=mutations[mutations['Impact']=='OFF']
        if homozygous_loss_only:
            loss = loss[loss['Zygosity']=='Homozygous']
        impacted = gain.append(loss)
        impacted = impacted.groupby(impacted.index).first()['Impact']
        try:    
            df
            df = pd.concat([df, pd.Series(impacted, name=cellline_name)], axis=1)
        except:
            print('done')
            df = pd.DataFrame(pd.Series(impacted, name=cellline_name))
    return df.transpose()


def transition_rates_calculation(df_rna_count_matrix, dict_nodes_genes=None,\
                                 max_percent = 0.95, diptest_max=0.05, bi_min=1.5, kurtosis_max=1,\
                                 amplification_factor = 100, return_initial_states=True):
    """ This function defines personalized transition rates values for Boolean model simulations with MaBoSS according to
        the workflow proposed in the publication of J. Béal et al. (2018), 
        'Personalization of Logical Models With Multi-Omics Data Allows Clinical Stratification of Patients' 
   
   ______ 
   inputs
    - df_rna_count_matrix (DataFrame): dataframe of RNA expression, where the samples are the lines and the genes are the columns.
    - dict_nodes_genes (dict):  dictionary, where the keys are the nodes of a model and the values the gene names corresponding to this nodes
    
    _______
    outputs
    - transition_rates_up (DataFrame): the transition rates of the activation of the genes for each sample. 
    They are sufficients to personalize the model but if the  transition_rates_down are needed
    the can be obtained by 1/tranisition_rates_up
    """
    df = df_rna_count_matrix.copy()
    if type(dict_nodes_genes)==dict:
        #filter the genes that are in the model
        dict_nodes2 = {node: dict_nodes_genes[node] for node in dict_nodes_genes if dict_nodes_genes[node] in df.columns}
        df = df.loc[:,dict_nodes2.values()]
        df = df.loc[:,~df.columns.duplicated()]


        #rename df with model nodes directly ?
        #inv_map = {v: k for k, v in dict_nodes2.items()}
        #df = df.rename(columns = inv_map)
        
    else:
        print("the dictionary providing the genes corresponding to the nodes of a model has not been defined")
    
    #filter the genes on the percent of cells their values are 0.
    df_filt = discard_percent_0s(df, max_percent=max_percent)
    
    #get the index of the genes that seems to be binarizable
    ind_bin = binarizable(df_filt, diptest_max=diptest_max, bi_min=bi_min, kurtosis_max=kurtosis_max)
    #normalzie these genes and store the result in bin_df
    print("binarization of {} genes started.".format(len(ind_bin)))
    bin_df = binarize(df_filt[ind_bin])
    print("binarization of {} genes done.\n".format(len(ind_bin)))

    #group all the other genes in another df
    no_bin_ind = [ind for ind in df_filt.columns if ind not in ind_bin]
    df_no_bin = df_filt[no_bin_ind]
    
    #get the index of the genes that seems to follow a 0 inflated distribution
    inflated0_ind = inflated0_test(df_no_bin)
    #normalize their values and store it in inflated0_df
    inflated0_df = inflated0_norm(df_no_bin[inflated0_ind])

    #get all the other indexes 
    univar_ind = [ind for ind in df_no_bin.columns if ind not in inflated0_ind]
    #normalize their expressions
    univar_df = uni_norm(df_no_bin[univar_ind])

    print("normalization of {} genes done.".format(len(no_bin_ind)))
    
    total_binnorm_df = pd.concat([bin_df, univar_df, inflated0_df], axis=1)
    #return the binarized/normalized values return total_binnorm_df

    transitions_up = amplification_factor**(2*(np.array(total_binnorm_df)-0.5))
    df_tr_up = pd.DataFrame(transitions_up, index=total_binnorm_df.index, columns=total_binnorm_df.columns)
    if return_initial_states == True:
        return df_tr_up, total_binnorm_df
    else:
        return df_tr_up

######################################################################
##                    MaBoSS simulation object                      ##
######################################################################

# Define the personalize_model function
def personalize_model(cell, cellname, model, dict_mut, dict_tr, dir, sample_count, previous):
    personalized_model = model.copy()
    nodes = model.network.names
    for node in nodes:
        # Assign node mutation
        if node in dict_mut:
            gene = dict_mut[node]
            if gene in cell['mutations'].keys():
                val_mut = cell['mutations'][gene]
                personalized_model.mutate(node, val_mut)
        # Assign transition rates
        if node in dict_tr:
            gene = dict_tr[node]
            if gene in cell['transition_rates_up'].keys():
                personalized_model.param['$u_'+node] = cell['transition_rates_up'][gene]
                personalized_model.param['$d_'+node] = 1/cell['transition_rates_up'][gene]
        # Assign initial states
        if len(cell['initial_states'])!=0:
            for node in cell['initial_states']:
                personalized_model.network.set_istate(node, cell['initial_states'][node], warnings= False)
                personalized_model.update_parameters(sample_count = sample_count)

    # Write the personal model in the sc_tmp directory
    personalized_model_path = os.path.join(dir,cellname)
    with open(personalized_model_path+'.bnd', 'w') as bnd_file:
        personalized_model.print_bnd(bnd_file)
    with open(personalized_model_path+'.cfg', 'w') as cfg_file:
        personalized_model.print_cfg(cfg_file)
    del personalized_model
    # For multiprocessing
    previous.append(cellname)

# Define the run simulation function
def run_cellmodel_simulation(cell, cellname, dir, previous):
    # Run and save the simulation results
    model_res = cell.run()
    model_res.save(dir+'/'+cellname+'_res')
    # For multiprocessing
    previous.append(cellname)

class CellEnsemble:
    def __init__(self, df_mutations, df_transition_rates, df_init_states, dict_gene_nodes, dict_strict_gene_nodes, project_name, df_groups):
        self.cells = {}
        self.mutations = df_mutations
        self.tr_rates = df_transition_rates
        self.init_cond = df_init_states
        self.dict_mut = dict_gene_nodes
        self.dict_tr = dict_strict_gene_nodes
        self.groups = df_groups
        # For project name -> If not defined - use date as project name
        if type(project_name) is str:
            self.name = project_name
        else:
            now = datetime.now()
            now = now.strftime("%Y_%m_%d_%Hh%M")
            self.name = print(now)

        # Transition rates - create a dictionary object for transition rates
        dic_tr = df_transition_rates.transpose().to_dict()
        
        # Initial conditions- create a dictionary object for initial conditions
        dic_init_states = df_init_states.transpose().to_dict()
        f_init_states = {}

        # For loop to create init cond for each cell
        for cell in dic_tr.keys():
            f_init_states[cell]={}
            for gene in dic_init_states[cell]:
                for node in dict_strict_gene_nodes:
                    if gene in dict_strict_gene_nodes[node]:
                        f_init_states[cell][node] = {0:1-dic_init_states[cell][gene], 1:dic_init_states[cell][gene]}                
        dic_init_states=f_init_states

        # For loop to create dictionary of parameters for each cells
        for cell in tqdm(dic_tr.keys()):

            # Get cellline name for each cell
            cellline = self.groups[cell]
            
            # Create mutation profile
            dic_mut = {k:v for k,v in df_mutations[cellline].to_dict().items() if v in ('ON','OFF')}
            
            # Create the dictionary of cells
            self.cells[cell] = {"name":cell,
                                "mutations":dic_mut,
                                "transition_rates_up":dic_tr[cell],
                                "dict_gene_nodes":dict_gene_nodes,
                                'initial_states':dic_init_states[cell],
                                'dict_strict_gene_nodes':dict_strict_gene_nodes
                                }

    def __repr__(self):
        return 'conditions : ' + str(list(self.model.keys())) + '\ngroups : ' + str(self.groups.cat.categories)

    def personalize_cellmodel(self, model, num_processes, reload_model = True, simulation_sample=1000,path=None):
        # Create base model condition
        self.model = {'base':{}}

        # Arg num_processes
        if num_processes == None:
            num_processes = 1
        elif type(num_processes) == int:
            num_processes = num_processes

        # Multiprocessing arg
        manager=mp.Manager()
        previous=manager.list()
        processes=[]

        # Set parameters
        self.core_model = model.copy()
        cells = list(self.cells.keys())

        # Create the sc_tmp directory
        condition_dir = None
        if path is None:
            base_folder = os.getcwd()
            if not os.path.exists(base_folder+'/sc_tmp/'):
                os.makedirs(base_folder+'/sc_tmp/')
                # Create the project directory
            project_folder = base_folder + '/sc_tmp/'+self.name
            if not os.path.exists(project_folder):
                os.makedirs(project_folder)
            else:
                now = datetime.now()
                now = now.strftime("%Y_%m_%d_%Hh%M")
                project_folder = project_folder + '_' + now
                os.makedirs(project_folder)
            self.project_dir = project_folder + '/'
            # Create base condition folder
            condition_dir = self.project_dir + 'base/'
            os.makedirs(condition_dir)
        else:
            if not os.path.exists(path):
                os.makedirs(path)
            self.project_dir = path
            condition_dir = path
            
        # Set the sample_count parameter
        if type(simulation_sample) == int:
            simulation_sample = simulation_sample


        # For loop for every group in the dataset
        print('Parameterizing the model')
        for group in tqdm(self.groups.cat.categories):
            
            # Create group directory
            group_dir = os.path.join(condition_dir,group)
            os.makedirs(group_dir, exist_ok=True)
            #print('Parameterizing model from group : ' + group + '\n')
            cells = self.groups[self.groups == group].index

            # Assign the sample_count for each cell group
            group_sample_count = len(cells)*simulation_sample
            
            # For loop for every cells in the group
            for i in range(len(cells)):
                cell = cells[i]
                while len(previous)<i-(num_processes-1):
                    time.sleep(1)
                p = mp.Process(target = personalize_model, 
                               args = (self.cells[cell], cell, self.core_model, self.dict_mut, self.dict_tr, group_dir, group_sample_count,previous))
                p.start()
                processes.append(p)
            for process in processes:
                process.join()
        
        # Reload the model
        if reload_model == True :
            print('loading model results...')
            for cell in tqdm(cells):
                self.model['base'][cell] = maboss.load(bnd_filename=group_dir+'/'+cell+'.bnd',
                                                       cfg_filenames=group_dir+'/'+cell+'.cfg')
                
    def load_personalized_model(self, dir):
        # Create base model condition
        self.model = {'base':{}}

        # Add the folder which 
        self.project_dir = dir+'/'
        cells = list(self.cells.keys())
        
        # Reload the model
        print('loading model results...')
        for cell in tqdm(cells):
            self.model['base'][cell] = maboss.load(bnd_filename=self.project_dir+cell+'.bnd',
                                                  cfg_filenames=self.project_dir+cell+'.cfg')

    def add_conditions(self, nodes, condition, condition_name, num_processes, simulation_sample=1000):
        # Arg num_processes
        if num_processes == None:
            num_processes = 1
        elif type(num_processes) == int:
            num_processes = num_processes

        # Multiprocessing arg
        manager=mp.Manager()
        previous=manager.list()
        processes=[]

        # Create base condition folder
        condition_dir = self.project_dir + condition_name
        os.makedirs(condition_dir)

        # Mutated the core model
        self.mutated_model = maboss.copy_and_mutate(self.core_model, nodes=nodes, mut = condition)

        # For loop for every group in the dataset
        print('Parameterizing the model')
        for group in self.groups.cat.categories:
            # Create group directories
            group_dir = condition_dir+'/'+group
            os.makedirs(group_dir)
            print('Parameterizing model from group : ' + group + '\n')
            cells = self.groups[self.groups == group].index

            # Assign the sample_count for each cell group
            group_sample_count = len(cells)*simulation_sample

            # For loop for every cells in the group
            for i in tqdm(range(len(cells))):
                cell = cells[i]
                while len(previous)<i-(num_processes-1):
                    time.sleep(1)
                p = mp.Process(target = personalize_model, 
                               args = (self.cells[cell], cell, self.mutated_model, self.dict_mut, self.dict_tr, group_dir, group_sample_count,previous))
                p.start()
                processes.append(p)
            for process in processes:
                process.join()

    def run_simulation(self, conditions='all', redo = False, num_processes = None):
        # Arg conditions
        if conditions == 'all':
            conditions = list(self.model.keys())
        elif type(conditions) not in [list, tuple]:
            conditions = [conditions]
        
        # Arg num_processes
        if num_processes == None:
            num_processes = 1
        elif type(num_processes) == int:
            num_processes = num_processes

        # Multiprocessing arg
        manager=mp.Manager()
        previous=manager.list()
        processes=[]
        
        # Create temporary results folder
        self.result_dir = self.project_dir + 'tmp_res/'

        # Get cells list
        cells = list(self.cells.keys())

        # For loop to run simulation for each condtion
        print('Compute simulations')
        for condition in conditions:
            print('\tComputing simulations for condition : ' + condition)
            condition_dir = self.result_dir + condition
            if not os.path.exists(condition_dir):
                os.makedirs(condition_dir)

            # For loop to run simulation for each cell
            for i in tqdm(range(len(cells))):
                cell = list(self.model[condition].keys())[i]
                if redo==False:
                    try:
                        self.model[condition][cell].result
                        print('the simulation {}|{} will not be re-computed'.format(cell, condition))
                        continue
                    except:
                        pass
                while len(previous)<i-(num_processes-1):
                    time.sleep(1)
                p = mp.Process(target = run_cellmodel_simulation, 
                               args = (self.model[condition][cell], cell, condition_dir, previous))
                p.start()
                processes.append(p)
            for process in processes:
                process.join()

    # [ ] - Not really necessary anymore since we have save modes in the sc_tmp folder
    def save_models(self, conditions='all'):
        # Create and define base directory for models
        ## Create scMODELs folder
        base_folder = os.getcwd()
        if not os.path.exists(base_folder+'/scMODELS'):
            os.makedirs(base_folder+'/scMODELS')
        ## Create project folder
        project_folder = base_folder+'/scMODELS/'+self.name
        if not os.path.exists(project_folder):
            os.makedirs(project_folder)
        else:
            now = datetime.now()
            now = now.strftime("%Y_%m_%d_%Hh%M")
            project_folder = project_folder + '_' + now
            os.makedirs(project_folder)
        self.project_dir = project_folder

        # Arg conditions
        if conditions == 'all':
            conditions = list(self.model.keys())
        elif type(conditions) not in [list, tuple]:
            conditions = [conditions]

        # For loop to run each conditions
        for condition in conditions:
            condition_dir = project_folder+'/'+condition
            os.makedirs(condition_dir)
            print('Saving models from condition:' + condition + '\n')
            # For loop for each cellline/samples
            for group in self.groups.cat.categories:
                group_dir = condition_dir+'/'+group
                os.makedirs(group_dir)
                cells = self.groups[self.groups==group].index
                print('\t Saving models from group : ' + group)
                for cell in tqdm(cells):
                    path_bnd = group_dir + '/' + cell + '.bnd'
                    path_cfg = group_dir + '/' + cell + '.cfg'
                    self.model[condition][cell].print_bnd(open(path_bnd,'w'))
                    self.model[condition][cell].print_cfg(open(path_cfg,'w'))

    def get_nodeprob_matrix(self, cells = 'all', conditions = 'all', fill_na = True):
        # Arg cell
        if cells == 'all':
            cells = list(self.results['base'].keys())
        elif type(cells) in [list, tuple]:
            cells = [cells]
        
        # Arg conditions
        if conditions == 'all':
            conditions = list(self.results.keys())
        elif type(conditions) not in [list, tuple]:
            conditions = [conditions]

        # for loop to obtain node_prob traj
        probtraj = {}
        for condition in conditions:
            probtraj_mtx=pd.DataFrame()
            for cell in cells:
                try:
                    probtraj_cell=self.results[condition][cell].get_last_nodes_probtraj()
                    probtraj_mtx = pd.concat([probtraj_mtx,probtraj_cell], ignore_index=True)
                except:
                    print('the simulation {}|{} seems not having been computed'.format(cells, condition))
            probtraj_mtx.index = cells
            if fill_na == True:
                probtraj_mtx = probtraj_mtx.fillna(0)
            probtraj[condition] = probtraj_mtx
        return probtraj
    
    def get_stateprob_matrix(self, cells = 'all', conditions = 'all', fill_na = True):
        # Arg cell
        if cells == 'all':
            cells = list(self.results['base'].keys())
        elif type(cells) in [list, tuple]:
            cells = [cells]
        
        # Arg conditions
        if conditions == 'all':
            conditions = list(self.results.keys())
        elif type(conditions) not in [list, tuple]:
            conditions = [conditions]
        
        # for loop to obtain state_prob traj
        probtraj = {}
        for condition in conditions:
            probtraj_mtx=pd.DataFrame()
            for cell in cells:
                try:
                    probtraj_cell=self.results[condition][cell].get_last_states_probtraj()
                    probtraj_mtx = pd.concat([probtraj_mtx,probtraj_cell], ignore_index=True)
                except:
                    print('the simulation {}|{} seems not having been computed'.format(cells, condition))
            probtraj_mtx.index = cells
            if fill_na == True:
                probtraj_mtx = probtraj_mtx.fillna(0)
            probtraj[condition] = probtraj_mtx
        return probtraj
#__________________________________________________________________________________________________________#
              
                
