import os
from tqdm import tqdm

class graph():
    def __init__(self,file, delimeter= ' '):

        self.rawgraph = file
        self.delimeter = delimeter
        self.num_nodes=0
        self.gene2id = {}

        # Set of unique genes in the graph
        self.genes = set()
        self.construct_graph()


    def construct_graph(self):
        self.graph= {}
        self.num_edges = 0
        with open(self.rawgraph, 'r') as f:
            edges = f.readlines()
            edges = [s.rstrip() for s in edges]

            for e in tqdm(edges):
                self.num_edges += 1
                src, dst,_ = e.split(self.delimeter)
                self.genes.update([src, dst])
                #src_, dst_ = self.getIDs(src, dst)
                if src in list(self.graph.keys()):
                    self.graph[src].append(dst)
                else:
                    self.graph[src]= [dst]

                if dst in list(self.graph.keys()):
                    self.graph[dst].append(src)
                else:
                    self.graph[dst]= [src]
                    
    def getIDs(self,src, dst):
        if src in list(self.gene2id.keys()):
            src = self.gene2id[src]
        else:
            self.gene2id[src] = self.num_nodes
            src = self.num_nodes
            self.num_nodes += 1

        if dst in list(self.gene2id.keys()):
            dst = self.gene2id[dst]
        else:
            self.gene2id[dst] = self.num_nodes
            dst = self.num_nodes
            self.num_nodes += 1
        return src,dst
    def getNodes(self):
         return list(self.genes)
    def summary(self):
        print('There are {} nodes, {} edges'.format(self.num_nodes, self.num_edges))

if __name__ == '__main__':

    g = graph('hintedge_new.txt')

    g.summary()

    print('Genes {}, Nodes {}'.format(len(g.genes),len(g.graph.keys())))
    print('{} \n{}'.format(list(g.genes)[:10], list(g.graph.keys())[:10]))
