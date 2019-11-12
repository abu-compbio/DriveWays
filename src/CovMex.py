import os
import numpy as np

class CovMex():
    def __init__(self, file):
        # print('New CovMex Created')
        self.gene2patients = {}
        self.file  = file
        self.read_data(self.file)

    def read_data(self,file):
        not_accepted = []
        #print('number of keys in graph: {}'.format(len(self.graph.gene2id.keys())))
        with open(file, 'r') as f:
            lines = f.readlines()
            lines = [s.rstrip() for s in lines]
            for line in lines:
                line = line.split(' ')
                # get the ID of the node from the graph object
                try:
                    g = line[0]#self.graph.gene2id[line[0]]
                except:
                    not_accepted.append(line[0])
                    continue
                # The gene already seen
                try:
                    # instead of doing if g in self.gene2patients.keys(),
                    # which may take O(n). I try to see if the key is not in the dict,
                    # if that's the case, there will be an error,
                    if isinstance(self.gene2patients[g], list):
                        raise Exception('There is an error in the file, genes appear in 2 different places')

                except:
                    self.gene2patients[g] = line[1:]
            self.num_patients = len(set([s for x in self.gene2patients.values() for s in x ]))
            #print('not accepted: ', len(not_accepted), len(lines))

    def getmutex(self,nodes):
        # print('type:', type(nodes))
        s= 0
        if isinstance(nodes, list):
            ps = []
            for node in nodes:
                try:
                    ps.extend(self.gene2patients[node])
                except:
                    continue
                    raise 'gene name {} is not in the gene_patient file'.format(node)

            # print(len(set(ps)),len(ps) )
            try:
                if len(ps) == 0: return 0.0 
                s = len(set(ps))/len(ps)
            except:
                print(nodes, len(set(ps)),len(ps))
            return s

        elif isinstance(nodes, str):
            # mutex of a single gene is 1
            return 1.0 if nodes in self.gene2patients.keys() else 0.0
            # return 1.0
        else:
            raise 'argument is not of type list or str'

    def getcov(self, nodes):
        # print('type:', type(nodes))
        if isinstance(nodes, list):

            ps = []
            for node in nodes:
                try:
                    ps.extend(self.gene2patients[node])
                except:
                    continue
                    raise Exception('gene name {} is not in the gene_patient file'.format(node))

            # print(len(set(ps)),len(ps) )
            return len(set(ps))/self.num_patients

        elif isinstance(nodes, str):

            return len(self.gene2patients[nodes])/self.num_patients if nodes in self.gene2patients.keys() else 0.0
        else:
            raise Exception('argument is not of type list or str')


    def summary(self):
        print('genes: ',len(self.gene2patients.keys()) )
        print('patients: ',len(set([s for x in self.gene2patients.values() for s in x ])))
        #print(' There are {} genes, and {} patients in {}'.format(len(self.gene2patients.keys()), len(set(self.gene2patients.values())),self.file ))


if __name__ == '__main__':
    from graph import Graph
    G = Graph('../data/intActedge_threshold_35.txt')
    m = CovMex('../data/gene_patient_update.txt')
    m.summary()
    # we know the cov and mutex of this modul from the original ClusterOne implementation
    # cov_ = 0.517685
    # mutex_ = 0.855928
    # module1 = ['TP53', 'CDKN2A' ,'PTEN', 'TP53INP1']
    # print(module1, '\n')
    # print('ClusterOne original implementation results\n COV: {}, MUTEX: {}'.format(cov_,mutex_))
    # print('COV: ',m.getcov(module1))
    # print('MUTEX: ',m.getmutex(module1))
