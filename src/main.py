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


class MexCoGrowth():

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


    # @profile
    def run(self):
        count_refuse = 0
        used_seeds = {}
        # total_genes= 0

        # This is added for a special MMR test on 14/05/2019
        total_genes_ = 0

        for seed in tqdm(self.seeds):
            try:
                # if a node is already in a cluster,it will not be used as a seed.
                _ = self.usednodes[seed]
                continue
            except:
                pass
            used_seeds[seed] = 1
            process = GrowthProcess(self.G,self.covmex,seed, self.usednodes, self.degrees, used_seeds, self.modules, self.strategy)
            module = process.grow()
            if module.size >= 3:
                self.modules.append(module)
                total_genes_ += len(module)
                # print('New Module: {}'.format(module.genes))
                total_genes = self.update_usednodes(module)

                self.update_degrees(module)
            else:
                # add the seed to used nodes, so it does not get used again later
                count_refuse +=1
                # print('Module: {} has size less than 3'.format(module.genes))

            if total_genes_ >= self.MAX_GENES:
                break
        # we pass the responsibility of generating top N genes file to the strategy
        # as the output of overlapping and no-overlapping is not the same
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


    graph_file = '../data/intActedge_threshold_35.txt'
    gene_patient_file = '../data/gene_patient_update.txt'
    seed_file = '../data/intAct_seeds.txt'
    condition_type = 'cov>cov/mutex'
    score_function_type = 'cov*mutex'


    #results_range = list(range(100,1200, 100))+[1173]+list(range(1200,1800, 100))+[1771]+list(range(1800,3400, 100))+[3368]
    # results_range = [720, 957, 418 ]

    results_range = [1771,3368,1173 ]



    # Randomization experiment
    # for i in range(100):
    #     print(f'{i+1}/100')
    #     components_folder = f'{date}/{i}/'
    #     seed_file = f'../data/random_seeds/intAct_seeds_threshold_35_{i}.txt'
    #     mkdir(components_folder)
    #
    #     algo = MexCoGrowth(
    #             graph_file,gene_patient_file,
    #             seed_file, condition_type, score_function_type, components_folder,
    #             contraction_allowed = True,
    #             MAX_GENES =  3368,t_threshold=1.4, d_threshold=4,results_range = results_range)
    #
    #     algo.run()

    # threshold paramas study
    # ds = [ 0,  1, 2, 3, 4,  5,6,  7,  8,9, 10]
    #ts = [0.8 , 0.92, 1.04, 1.16, 1.28, 1.4 , 1.52, 1.64, 1.76, 1.88, 2.  ]
    #ts = [0.32, 0.44, 0.56] #[ 0.32, 0.44, 0.56]#,
    #ts = ts[::-1]
    # for d in ds:
    # for t in ts:
    print('t {}/d {}'.format(t,d))
    components = '../out/components/'
    date = components+'24_06_2019'

    components_folder = '{}/d{}_t{}/'.format(date,d,t )
    mkdir(components_folder)
    algo = MexCoGrowth(
            graph_file,gene_patient_file,
            seed_file, condition_type, score_function_type, components_folder,
            contraction_allowed = True,
            MAX_GENES =  3368,t_threshold=t, d_threshold=d,results_range = results_range)

    algo.run()



    # # Cancer specific study

    # graph_file = '../data/intActedge_threshold_35_BRCA.txt'
    # gene_patient_file = '../data/gene_patient_BRCA.txt'
    # seed_file = '../data/intAct_seeds_BRCA.txt'
    # condition_type = 'cov>cov/mutex'
    # score_function_type = 'cov*mutex'
    # components_folder = f'{date}/BRCA/'
    # mkdir(components_folder)
    #
    # algo = MexCoGrowth(
    #         graph_file,gene_patient_file,
    #         seed_file, condition_type, score_function_type, components_folder,
    #         contraction_allowed = True,
    #         MAX_GENES =  1000,t_threshold=t, d_threshold=d,results_range = results_range)
    #
    # algo.run()


if __name__ == '__main__':

    ds = [ 1, 2, 3, 4,  5]
    ts = [0.32, 0.44, 0.56, 0.68]
    ts = ts[::-1]
    # # p = Pool(5)
    for t in ts:
        for d in ds:
            main(t,d)

    # main()





    # i_ = 0
    # seeds = ['seeds_cov','seeds']
    # seeds_name = ['cov_seed', 'cov*mutex_seed']
    # conditions_ = ['cov*mutex'] #['cov/mutex', 'cov>cov/mutex'] #, 'cov*mutex'
    # conditions_names = ['cov*mutex'] #['cov_over_mutex', 'cov_2'] #, 'cov*mutex'
    # scores_ =['cov', 'cov*mutex']
    # for seed_, seed_name in zip(seeds, seeds_name):
    #     seed_file = '../data/{}.txt'.format(seed_)
    #     for cond_, cond_name in zip(conditions_, conditions_names):
    #         for score_ in scores_:
    #             # if i_ <4:
    #             #     i_ +=1
    #             #     continue
    #             components_folder = '{}/{}_C_{}_S_{}/'.format(date,seed_name,cond_name,score_ )
    #             mkdir(components_folder)
    #             print('\n', '{}_C_{}_S_{}/'.format(seed_name,cond_name,score_ ))
    #             algo = ClusterOne(
    #                     graph_file,gene_patient_file,
    #                     seed_file, cond_, score_, components_folder,
    #                     contraction_allowed = True,
    #                     MAX_GENES = 2500)
    #
    #             algo.run()
