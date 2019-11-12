# This class takes a seed and usednodes,
#and follow a startegy to grow a module

from strategy.condition import Condition
from strategy.NoOverlap import NoOverlapStrategy
from strategy.functions import *
from strategy import *
from  config import *
from module import Module
from CovMex import *
from graph import Graph

class GrowthProcess():
    def __init__(
            self,graph,cov_mutx,
            seed, usednodes, degrees, used_as_seed, prev_modules,strategy):

        self._G = graph
        self.seed= seed
        self.usednodes = usednodes
        self.degrees = degrees
        self.prev_modules = prev_modules
        self.strategy = strategy
        self.module = Module(cov_mutx, seed)
        self.used_as_seed = used_as_seed

    @property
    def graph(self):
        return self._G
    @graph.setter
    def graph(self, G):
        self._G = G

    def grow(self):
        should_stop = False
        while( not should_stop):
            current_module = self.module.copy()
            # print('Current module: ', current_module.genes)
            externel_nodes = self._G.getBoundary_nodes(current_module.genes)
            articulation_points = self._G.get_articulation_points(current_module.genes)

            # print('externel_nodes: ', externel_nodes)

            nodes, should_add, should_stop = self.strategy(current_module,externel_nodes, articulation_points,
                                    self.usednodes, self.degrees, self._G.G, self.used_as_seed, self.prev_modules)

            if should_stop:
                # print('module: {}, current cov {}, and mutex {}'.format(self.module.genes, self.module.cov, self.module.mutex))
                break
            if should_add:
                # print(nodes, ' are added. ')
                self.module.add_gene(nodes)
                # print('current cov {}, and mutex {}'.format(self.module.cov, self.module.mutex))
            else:
                self.module.remove_gene(nodes)


        return self.module

if __name__ == '__main__':
    G = Graph('../data/hintedge_new.txt')
    covmex = CovMex('../data/gene_patient_update.txt')
    module1 = []#['TP53', 'CDKN2A' ,'PTEN', 'TP53INP1']
    usednodes = {k:1 for k in module1}
    # print(covmex.getcov('RNF14'))
    nodes = G.nodes
    # m = Module(covmex, nodes[0])
    condition = Condition(threshold=0.9)
    score = cov
    strategy = NoOverlapStrategy(condition, score)
    process = GrowthProcess(G,covmex,'TP53INP1', usednodes, strategy)
    result = process.grow()

    # module1 = ['TP53', 'CDKN2A' ,'PTEN', 'TP53INP1']
    # c1 = covmex.getcov(module1)
    # m1 = covmex.getmutex(module1)
    # print(module1,': COV: ',c1, ' MUTEX: ',m1, ', COV/MUTEX: ', c1/m1)
    # c2 = result.cov
    # m2 = result.mutex
    # print(result.genes,': COV: {}, MUTX: {}, COV/MUTEX: '.format(c2,m2,c2/m2))
    # print('cov/mutex *t: {} > covinit/ mutexinit: {}'.format(c1/m1*0.9, c2/m2) )

    print(result)
