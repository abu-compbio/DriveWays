# TODO:
#    - Understand the mapping between gene name and ID in cluster_one
#    - Generate a file where is line contains the seed
# the seed are soted in descending order based on avg_cov + mutex
# the graph is defined as an adjacency list
import os
import numpy as np
from tqdm import tqdm
from MexCov import *
from graph import *


class SeedGenerater():

    def __init__(self, ppi_file, gene_patient_file, sep=' '):
        self.ppi_file  = ppi_file
        self.gene_patient_file = gene_patient_file
        self.graph = graph(self.ppi_file,sep)
        self.MexCov = MexCov(self.gene_patient_file, self.graph)
        self.seeds = self.generate()



    def generate(self):
        nodes = list(self.graph.getNodes())
        print('number of nodes: ', len(nodes))
        scores = [self.MexCov.get_score(n) for n in nodes]
        sorted_scores_ids = np.argsort(scores)[::-1]
        sorted_scores = np.array(scores)[sorted_scores_ids]
        # print(sorted_scores[:20])
        sorted_seeds = np.array(nodes)[sorted_scores_ids]

        return sorted_seeds

    def random_sort(self, n):
        nodes = list(self.graph.getNodes())
        seeds_lists = []
        for i in range(n):
            seeds = np.array(nodes)
            np.random.shuffle(seeds) # inplace shuffle
            seeds_lists.append(seeds)

        return seeds_lists

    def write_seeds(self, seed_file):
        with open(seed_file, 'w') as f:
            for s in self.seeds:
                f.write(s+'\n')



if __name__ == '__main__':
    generator = SeedGenerater('../../data/intActedge_threshold_35_BRCA.txt', '../../data/gene_patient_BRCA.txt', '\t')
    seeds = generator.seeds
    print(seeds)
    generator.write_seeds('intAct_seeds_BRCA.txt')
    # randomized_lists = generator.random_sort(100)
    # for i in range(100):
    #     with open(f'../../data/random_seeds/intAct_seeds_threshold_35_{i}.txt', 'w') as f:
    #         for s in randomized_lists[i]:
    #             f.write(s+'\n')
