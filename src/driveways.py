import os
import sys
import time
import numpy as np
from tqdm import tqdm
import glob
from strategy.condition import Condition
from strategy.NoOverlap import NoOverlapStrategy
from strategy.Overlap import OverlapStrategy
from strategy.functions import *
from strategy import *
from  config import *
from module import Module
from CovMex import *
from graph import Graph
from growthprocess import GrowthProcess


from multiprocessing import Pool


class OLDRIM():

    def __init__(
            self,graph_file,gene_patient_file,
            seed_file, condition_type, score_function_type,
            components_folder,Process_seed =False,
            contraction_allowed = True,MAX_GENES = 2500,t_threshold=0.9, d_threshold=1, results_range=None):

        if not isinstance(condition_type, str) or not isinstance(score_function_type, str) :
            raise Exception('The condition_type and/or score_function_type should be of type str')
        if not condition_type in fncs.keys() and not score_function_type in fncs.keys() :
            raise Exception('The condition_type and score_function_type should be one of the following: {}'.format(list(functions.keys())))
        self.G = Graph(graph_file)
        self.covmex = CovMex(gene_patient_file)
        # consider removing seed during the grwoth process
        self.Process_seed = Process_seed
        # Consider removing genes as well, not just adding
        self.contraction_allowed = contraction_allowed
        self.usednodes = {}
        self.degrees = {}
        self.condition = Condition(condition_type, threshold=t_threshold)
        self.score_fnc = fncs[score_function_type]
        self.strategy = OverlapStrategy(self.condition, self.score_fnc, components_folder, self.Process_seed, self.contraction_allowed,d_threshold)
        self.seeds = [s.strip() for s in open(seed_file).readlines()]
        # process = GrowthProcess(self.G,self.covmex,'TP53INP1', self.usednodes, self.strategy)
        self.modules = []
        self.MAX_GENES = MAX_GENES
        self.results_range = results_range



    def run(self):
        count_refuse = 0
        used_seeds = {}
        total_genes = 0

        for seed in tqdm(self.seeds):
            try:
                # if a node is already in a module,it will not be used as a seed.
                _ = self.usednodes[seed]
                continue
            except:
                pass
            used_seeds[seed] = 1
            process = GrowthProcess(self.G,self.covmex,seed, self.usednodes,
                                    self.degrees, used_seeds,
                                    self.modules, self.strategy)

            module = process.grow()

            if module.size >= 3:
                self.modules.append(module)
                total_genes += len(module)
                self.update_degrees(module)
            else:
                # add the seed to used nodes, so it does not get used again later
                count_refuse +=1

            if total_genes >= self.MAX_GENES:
                break

        print('total_genes: ', total_genes)
        print('{} modules were detected. \n {} module were refused'.format(len(self.modules), count_refuse))
        self.strategy.post_process_modules(self.modules, self.MAX_GENES, self.results_range)
        return

    def update_usednodes(self,module):
        for g in module.genes:
            try:
                # len(self.modules)-1 is the index of the module to which the gene g belongs
                if self.usednodes[g]:
                    self.usednodes[g].append(len(self.modules)-1)
            except:
                self.usednodes[g] = [len(self.modules)-1]

        return len(list(self.usednodes.keys()))

    def update_degrees(self,module):

        for g in module.genes:
            try:
                # len(self.modules)-1 is the index of the module to which the gene g belongs
                if self.degrees[g]:
                    deg_ = len(set(list(self.G.G.adj[g])).intersection(set(module.genes)))
                    self.degrees[g].append(deg_)
            except:
                deg_ = len(set(list(self.G.G.adj[g])).intersection(set(module.genes)))
                self.degrees[g] = [deg_]

        return






def main(t =1.4, d=4):

    algo = OLDRIM(
            graph_file,gene_patient_file,
            seed_file, condition_type, score_function_type, components_folder,
            contraction_allowed = True, MAX_GENES =  3368,t_threshold=t,
            d_threshold=d,results_range = results_range)

    algo.run()




if __name__ == '__main__':

    t = 1.4
    d = 4

    main(t,d)
