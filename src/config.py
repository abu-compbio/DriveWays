from strategy.functions import *

graph_file = '../data/intAct_PPI.txt'
gene_patient_file = '../data/gene_patients.txt'
seed_file = '../data/intAct_seeds.txt'
condition_type = 'cov>cov/mutex'
score_function_type = 'cov*mutex'

#for every total_genes in results_range, the algorithm will output the identified modules to a file in out/ folder,
# under the name cc_n<N>_k3.txt, each line represents a module.
# The values here are the total_genes for KEGG, Reactome, and BioCarta, respectively as described in the paper.
results_range = [3368,1173,1771 ]
