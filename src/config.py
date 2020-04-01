from strategy.functions import *

graph_file = '../data/intAct_PPI.txt'
gene_patient_file = '../data/gene_patients_Pan-cancer.txt'
seed_file = '../data/intAct_seeds.txt'
condition_type = 'cov>cov/mutex'
score_function_type = 'cov*mutex'
# min size
k=3

#for every total_genes in results_range, the algorithm will output the identified modules to a file in out/ folder,
# under the name cc_n<N>_k3.txt, each line represents a module.
# The value here is the total_genes for KEGG as described in the paper.
results_range = [1173 ]


# GP Parameters
t_low,t_high = [0.8,1.2]
d_low,d_high = [2,5]
GP_n_calls = 11
seed = 1234
GP_verbose = True
GP_n_jobs = -1
