# This class takes a seed and usednodes,
#and follow a startegy to grow a module

from strategy.condition import Condition
from strategy.Overlap import OverlapStrategy
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

            externel_nodes = self._G.getBoundary_nodes(current_module.genes)
            articulation_points = self._G.get_articulation_points(current_module.genes)

            nodes, should_add, should_stop = self.strategy(current_module,externel_nodes, articulation_points,
                                    self.usednodes, self.degrees, self._G.G, self.used_as_seed, self.prev_modules)

            if should_stop:
                break
            if should_add:
                self.module.add_gene(nodes)
            else:
                self.module.remove_gene(nodes)

        return self.module

if __name__ == '__main__':
    PPI_file = '../data/intAct_PPI.txt'
    mutation_data_file = '../data/gene_patients.txt'
    G = Graph(PPI_file)
    covmex = CovMex(mutation_data_file, g)
    # This is a toy example to check that the growth process is properly implemented
    module1 = ['TP53', 'CDKN2A' ,'PTEN', 'TP53INP1']
    usednodes = {k:1 for k in module1}
    nodes = G.nodes

    condition = Condition(threshold=0.9)
    score = cov
    strategy = OverlapStrategy(condition, score)
    process = GrowthProcess(G,covmex,'TP53INP1', usednodes, strategy)
    result = process.grow()

    print(result)
