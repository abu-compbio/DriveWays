import os
import sys
import time
import numpy as np
from tqdm import tqdm
import glob
from strategy.condition import Condition
from strategy.Overlap import OverlapStrategy
from strategy.functions import *
from strategy import *
from  config import *
from module import Module
from CovMex import *
from graph import Graph
from growthprocess import GrowthProcess


from multiprocessing import Pool
import skopt
from skopt import gp_minimize, callbacks
from skopt.space import Real, Integer, Categorical
from skopt.utils import use_named_args
from  skopt.plots import plot_convergence


class DriveWays():

    def __init__(
            self,graph_file,gene_patient_file,
            seed_file, condition_type, score_function_type,
            components_folder,Process_seed =False,
            contraction_allowed = True,MAX_GENES = 2500,UNIQUE_GENES = None,
            t_threshold=0.9, d_threshold=1, results_range=None):

        if not isinstance(condition_type, str) or not isinstance(score_function_type, str) :
            raise Exception('The condition_type and/or score_function_type should be of type str')
        if not condition_type in fncs.keys() and not score_function_type in fncs.keys() :
            raise Exception('The condition_type and score_function_type should be one of the following: {}'.format(list(functions.keys())))
        if not MAX_GENES and not UNIQUE_GENES:
            raise Exception('Either MAX_GENES or UNIQUE_GENES should be set, the other should be None')

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
        self.UNIQUE_GENES = UNIQUE_GENES
        self.results_range = results_range


    # @profile
    def run(self):
        count_refuse = 0
        used_seeds = {}
        total_genes = 0
        unique_genes = 0
        seed_skipped = 0
        for seed in tqdm(self.seeds):
            if seed in self.usednodes.keys():
                seed_skipped += 1
                continue
            # Add the seed to used_seeds
            used_seeds[seed] = 1
            process = GrowthProcess(self.G,self.covmex,seed, self.usednodes, self.degrees, used_seeds, self.modules, self.strategy)
            module = process.grow()
            if module.size >= 3:
                self.modules.append(module)
                total_genes += len(module)
                # print('New Module: {}'.format(module.genes))
                unique_genes = self.update_usednodes(module)

                self.update_degrees(module)
            else:
                count_refuse +=1

            if self.MAX_GENES:
                if total_genes >= self.MAX_GENES:
                    break
            else:
                if unique_genes >= self.UNIQUE_GENES:
                    break

        post_modules =  self.strategy.post_process_modules(self.modules, self.MAX_GENES,self.UNIQUE_GENES, self.results_range)
        ODMSS = 0
        for pred_m in post_modules:
            ODMSS += self.covmex.getcov(pred_m)* self.covmex.getmutex(pred_m)
        return ODMSS

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


t = Real(low=t_low, high= t_high, prior='log-uniform', name='t_threshold')
d = Real(low=d_low, high = d_high, prior = 'log-uniform', name='d_threshold')
default_parameters = [1.0,3]

dimensions = [t,d]
@use_named_args(dimensions=dimensions)
def get_validation_score(t_threshold,d_threshold):

    #print(init_learning_rate)
    t = np.float32(t_threshold)
    d = np.float32(d_threshold)

    result_folder = '../out/components/GP/'
    mkdir(result_folder)
    components_folder = f'{result_folder}/d{d}_t{t}/'
    mkdir(components_folder)
    algo = DriveWays(
            graph_file,gene_patient_file,
            seed_file, condition_type, score_function_type, components_folder,
            contraction_allowed = True,
            MAX_GENES =  1771,UNIQUE_GENES = None, t_threshold=t, d_threshold=d,results_range = results_range)

    ODMSS = algo.run()



    return -ODMSS


def main():

    save_folder = '../out/components/GP/'
    checkpoint_callback = callbacks.CheckpointSaver(f'{save_folder}/search_result.pkl')

    search_result = gp_minimize(func=get_validation_score,
                            dimensions=dimensions,
                            acq_func='EI',  # Expected Improvement.
                            n_calls=GP_n_calls,
                            x0=default_parameters,
                            random_state=seed, verbose =GP_verbose,
                            callback = [checkpoint_callback],
                           n_jobs=GP_n_jobs)

    import shutil
    best_t,best_d = search_result.x
    # a trick to get the folder containing corresponding to the best result
    t,d = str(best_t), str(best_d)
    t,d = float(t[:5]), float(d[:5])
    try:
        src = glob.glob(f'../out/components/GP/d{d}*_t{t}*/')[0]
        dst = f'../out/components/DriveWays/'
    except:
        pritn('src ', src)
    shutil.copytree(src, dst)
    print(f'results are saved to /out/components/DriveWays/ folder')


if __name__ == '__main__':
    main()
