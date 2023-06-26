
import os
import numpy as np
import pandas as pd
import maboss

import pyprofile.pyPROFILE_mutation as mut
import pyprofile.pyPROFILE_v0_4_6 as profile

import argparse
 
parser = argparse.ArgumentParser(description="Building PhysiCell simulation for a particular ",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-r", "--rnaseq", help="RNASeq data file", required=True)
parser.add_argument("-m", "--mutations", help="Mutations data file", required=True)
parser.add_argument('-l', "--cell_lines", help="File with the list of cell lines to use")
parser.add_argument("-b", "--bnd_file", help="MaBoSS BND model file", required=True)
parser.add_argument("-c", "--cfg_file", help="MaBoSS CFG model file", required=True)
parser.add_argument("-d", "--dictionnary", help="Dictionary associating MaBoSS nodes to Genes ID", required=True)

args = parser.parse_args()

CL_RNA = pd.read_csv(args.rnaseq, sep = ',', index_col=1)
CL_RNA.index.name = None

CL_mut = pd.read_csv(args.mutations, sep = ',')
CL_mut['model_name'] = CL_mut['model_name'].astype('category')

if args.cell_lines is None:
    common_CL = list(set(CL_RNA.index) & set(CL_mut['model_name'].cat.categories))
else:    
    with open(args.cell_lines, 'r') as cl_file:
        common_CL = [cl.strip() for cl in cl_file.readlines()]

CL_RNA = CL_RNA[CL_RNA.index.isin(common_CL)]
CL_mut = CL_mut[CL_mut.model_name.isin(common_CL)]

bnd_file = args.bnd_file
cfg_WT = args.cfg_file

## Load MaBoSS model
sizek = maboss.load(bnd_file, cfg_WT)

fname = 'Models/'
# Dictionary file for the node transition rates and initial conditions
dict_strict_gene_nodes=pd.read_csv(args.dictionnary, index_col=None)
dict_strict_gene_nodes=dict(np.array(dict_strict_gene_nodes))

# Dictionary file for the node mutations
dict_gene_nodes=pd.read_csv(args.dictionnary, index_col=None)
dict_gene_nodes=dict(np.array(dict_gene_nodes))


mutation_data = mut.mutation_assignment(CL_mut, database_dir = 'data')
mutation_data = mut.mutation_fusion(mutation_data)
mutation_data.index.name = None
mutation_data.columns.name = None

# Replace Boolean mutation with ON and OFF
mutation_data = mutation_data.replace(0.0, 'OFF')
mutation_data = mutation_data.replace(1.0, 'ON')

# Calculate initial condition and transition rates
sc_rates = profile.transition_rates_calculation(CL_RNA, 
                                                 dict_nodes_genes=dict_strict_gene_nodes,
                                                 return_initial_states = True,
                                                 amplification_factor = 30)

cell_line_index = pd.Series(CL_RNA.index, index = CL_RNA.index).astype('category')
print(cell_line_index)
test_scPROFILE = profile.CellEnsemble(df_mutations= mutation_data,
                                df_transition_rates= sc_rates[0],
                                df_init_states=sc_rates[1],
                                dict_gene_nodes=dict_gene_nodes,
                                dict_strict_gene_nodes=dict_strict_gene_nodes,
                                project_name = 'CL_sizek',
                                df_groups=cell_line_index)

test_scPROFILE.personalize_cellmodel(sizek, 
                                    num_processes = 20,
                                    reload_model = False,
                                    path="models/personalized")




