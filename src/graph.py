import networkx as nx
from networkx.algorithms import boundary as bn
from CovMex import *

class Graph():
    def __init__(self, file, delimeter= ' '):

        self.rawgraph = file
        self.delimeter = delimeter
        self.construct_graph()


    def construct_graph(self):
        self.G= nx.Graph()
        self.num_edges = 0
        # mutexs= []
        with open(self.rawgraph, 'r') as f:
            edges = f.readlines()
            edges = [tuple(s.strip().split('\t')[:2]) for s in edges]
            # mutexs = [covmex.getmutex(list(edge)) for edge in edges]
            # print(edges[:3])
            self.G.add_edges_from(edges)
        self.G = self.G.to_undirected()
        # return mutexs,edges

    @property
    def size(self):
        return self.G.number_of_nodes()
    @property
    def NumberEdges(self):
        return self.G.number_of_edges()
    @property
    def nodes(self):
        return list(self.G.nodes)

    def getDegree(self,gene=None):
        if gene:
            return self.G.degree[gene]
        else:
            return self.G.degree


    def getEdges(self,val):
        if isinstance(val, list):
            return slef.G.edges(val)
        elif isinstance(val, str) :
            return self.G.edges(val)
        return list(self.G.edges)

    def getAdj(self,val):
        if isinstance(val, list):
            # dict_edges = self.G.edges(val)
            # adjs = set(dict_edges.values())
            return list(self.G.edges(val))
        elif isinstance(val, str) :
            return list(self.G[val])

        return list(self.G.edges)

    # Return a set of nodes that are adjecent to the nodes in nodes
    def getBoundary_nodes(self,nodes, out_nodes=None):
        return bn.node_boundary(self.G, nodes,out_nodes )
    def get_articulation_points(self, nodes):
        sub_G = self.G.subgraph(nodes)
        points = nx.articulation_points(sub_G)
        return list(points)

    def summary(self):
        print('There are {} nodes, {} edges'.format(self.size, self.NumberEdges))


if __name__ == '__main__':
    # print('nx version: ', nx.__version__)
    import numpy as np
    from scipy import stats
    G = Graph('../data/intActedge_threshold_35.txt')


    # m = CovMex('../data/gene_patient_update.txt')
    # mutexs, edges = G.construct_graph(m)
    # mutexs, edges = np.array(mutexs), np.array(edges)
    # print(stats.describe(mutexs))
    # print('Qunatiles: \n', np.quantile(mutexs, q=[0.05,0.1,0.2,0.25]))
    # print('edges mutex < 0.9: \n', len(mutexs[mutexs<0.95]))
    # s_edges = edges[mutexs<0.95]
    # print(s_edges)
    # tp53 = [1 for e in s_edges if e[0] == 'TP53' or e[1] == 'TP53']
    # print('TP53s: ', sum(tp53))

    # nodes = G.getNodes()
    # edges = G.getEdges(None)
    G.summary()
    # b = bn.node_boundary(G.G, nodes[:2])
    # print(b)
    # print(G.getAdj(nodes[0]))
    # print(G.getAdj(nodes[:5]))
    # print('Adj of {}: {}'.format(nodes[0], G.getAdj(nodes[0])))
    # print('Adj of {}: {}'.format(nodes[1], G.getAdj(nodes[1])))
    # print('Adj of {}: {}'.format(nodes[:2], G.getAdj(nodes[:2])))

    # a= ['COL4A2', 'HDLBP', 'BAG6', 'PSG3', 'TBX4', 'ATXN1', 'EGLN2', 'SPINT1', 'KEAP1', 'NR2F1', 'NR2F2', 'ZFYVE9', 'BCCIP', 'HTRA1']
    # c = ['UBE2E1', 'UBE2D4', 'UBE2D2', 'UBE2D3', 'UBE2E3', 'UBE2D1', 'AR', 'DACH1', 'UBE2V1', 'UBE2E2', 'UBE2U', 'UBE2W']
    #
    # s = set(a) | set(c)
    # print(len(s), len(b),len(b.intersection(s)))
