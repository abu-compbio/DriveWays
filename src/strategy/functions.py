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
