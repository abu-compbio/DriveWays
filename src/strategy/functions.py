# This is
import os
import numpy as np
import cProfile, pstats, io
import networkx as nx
from networkx.algorithms import matching
def cov_over_mutex(module, t=None):

    return np.float128(module.cov)/ np.float128(module.mutex) if t== None else (np.float128(module.cov)/ np.float128(module.mutex))*np.float128(t)

def cov(module, t=None):
    return np.float128(module.cov) if t== None else np.float128(module.cov)*np.float128(t)

def mutex(module, t=None):
    return np.float128(1) if t== None else np.float128(module.mutex)/np.float128(t)

def cov_mutx(module, t=None):
    return np.float128(module.cov)* np.float128(module.mutex) if t== None else (np.float128(module.cov)* np.float128(module.mutex))*np.float128(t)

def cov_cov_mutex(module, t=None):
    # This is similar to the buggy implementation in java, discovered on 08/02/2019

    return np.float128(module.cov)/ np.float128(module.mutex)  if t== None else np.float128(module.cov)/np.float128(t)


fncs = {
        'cov/mutex': cov_over_mutex,
        'cov': cov,
        'mutex':mutex,
        'cov*mutex': cov_mutx,
        'cov>cov/mutex': cov_cov_mutex
}


def mkdir(dir):
    if not os.path.exists(dir):
        os.mkdir(dir)




def predict_modules(file):
    with open(file) as f:
        lines = f.readlines()
    modules = [l.strip().split() for l in lines]
    modules_ = []
    mapping = {}
    index = 0
    for i,l in enumerate(modules):
        mapping[index] = '{}'.format(i)
        index +=1
        modules_.append(l)
    return modules_, mapping

def txt_parser2(file):

    with open(file) as f:
        lines = f.readlines()
        modules = [l.strip().split() for l in lines]

    return modules

def read_gene_assocition(file):
    association = {}
    with open(file) as f:
        for l in f.readlines():
            l = l.strip().split()
            association[l[0]] = l[1:] if len(l) > 1 else []
    return association

def get_mmr(edges,n_ref, n_pred):
    G= nx.Graph()
    G.add_weighted_edges_from(edges)
    G = G.to_undirected()

    max_match = matching.max_weight_matching(G, maxcardinality=True)
    sum_ = 0
    for e in max_match:
        sum_ += G[e[0]][e[1]]['weight']
    return sum_/min(n_ref, n_pred)


from scipy.stats import hypergeom
from statsmodels.stats import multitest


# Load IntAct genes, the number of genes in the network is used as the background
# value for q-value calculation
with open('../data/intact_genes.txt') as f:
    intAct_genes = [s.strip().split('\t')[-1] for s in f.readlines()]

#Load the GO standard association dictionary, this is used to calculated MMR with GO consistency score
level5_assocition = read_gene_assocition('../data/go_term_intact_association_level_5.txt')


def overlap(ref, pred):
    return len(set(ref).intersection(set(pred)))**2/(len(set(ref))*len(set(pred)))
def overlaps_func(ref_modules,pred_modules):
    scores = []
    values = []
    for i,ref_module in enumerate(ref_modules):
        matches = []
        for j, pred_module in  enumerate(pred_modules):
            scores.append(overlap(ref_modules[i], pred_modules[j]))
            #values.append(len(set(ref_modules[i]).intersection(set(pred_modules[j]))))
    return scores#, values

def qvalues(ref_modules,pred_modules):
    M = len(intAct_genes)
    pvals = []
    for i in range(len(ref_modules)):
        for j in range(len(pred_modules)):
            a = ref_modules[i]
            b = pred_modules[j]
            n = len(a)
            N = len(b)
            x = len(set(a).intersection(set(b)))
            pvals.append(hypergeom.sf(x-1, M, n, N))

    qvals = multitest.multipletests(pvals, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)

    return qvals[1]

def GO_Consistency_score(ref_modules,pred_modules):

    GO_cons_level5 = []

    for i in range(len(ref_modules)):#
        ref_GO_terms_level5 = set().union(*[set(level5_assocition[g]) for g in ref_modules[i] if g in level5_assocition])



        for j in range(len(pred_modules)):#
            pred_GO_terms_level5 = set().union(*[set(level5_assocition[g]) for g in pred_modules[j] if g in level5_assocition])
            if len(set(ref_GO_terms_level5)) == 0:
                w_level5 = 0
            else:
                w_level5 = len(set(ref_GO_terms_level5).intersection(set(pred_GO_terms_level5)))/len(set(ref_GO_terms_level5).union(set(pred_GO_terms_level5)))
            GO_cons_level5.append(w_level5)

    return GO_cons_level5
