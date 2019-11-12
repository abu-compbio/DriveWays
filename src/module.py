## This is a class that represent a module
# It will have functions, like calculating the coverage, mutex and so on
# it will have a function to caclulate the overlapp with other module

from CovMex import *

class Module():

    def __init__(self,cov_mutex, seed=None):
        '''
        genes (list):
        cov_mutex (COVMEX): object for calculating the coverage and mutex of a given module
        '''
        self._genes = []
        self._seed = None
        if isinstance(seed, list):
            # This condition appeared when we call copy function on an instance of this class
            self._genes = seed
        elif isinstance(seed, str):
            self._genes.append(seed)
            self._seed= seed

        self._seed = seed
        self.cov_mutex =cov_mutex

    def is_seed(self, node):
        return node == self._seed
    @property
    def size(self):
        return len(self._genes)

    @property
    def cov(self):
        return self._cov(self)

    def _cov(self,module = None):
        if not isinstance(module, Module):
            raise 'Module should be of instance module'
        # print('module.genes: ',module.genes)
        return self.cov_mutex.getcov(module.genes)

    @property
    def mutex(self):
        return self._mutex(self)

    def _mutex(self,module = None):
        if not isinstance(module, Module):
            raise Exception('Module should be of instance module')
        return self.cov_mutex.getmutex(module.genes)

    @property
    def genes(self):
        return self._genes[:]

    @genes.setter
    def genes(self, genes):
        self._genes = genes

    # dunder functions

    def __len__(self):

        return len(self.genes)

    def __str__(self):

        return ' '.join(self._genes)

    def __contains__(self, gene):

        return gene in self._genes

    def __getitem__(self, key):

        return self._genes[key]

    def copy(self):
        return Module(self.cov_mutex,self.genes[:])

    def add_gene(self,genes):
        if isinstance(genes, str):
            self._genes.append(genes)
        elif isinstance(genes, list):
            self._genes.extend(genes)
        else:
            raise Exception('genes should be either of type str or list')
        return self
        
    def remove_gene(self,genes):
        l = len(self.genes)
        if isinstance(genes, str):
            self._genes.remove(genes)
        elif isinstance(genes, list):
            self._genes = [s for s in self.genes if s not in genes]
        if l == len(self.genes):
            raise Exception('Requested genes to be removed are not in the module')
        return self

    def overlap(self, module, threshold=0.2):
        a = len(set(self.genes).intersection(set(module.genes)))/len(set(self.genes).union(set(module.genes)))
        return a < threshold






if __name__ == '__main__':
    covmex_ = CovMex('../data/gene_patient.txt')
    covmex_.summary()
    m = Module(covmex_)
    m.add_gene('geneA')
    m.add_gene(['geneB','geneC'])
    print(m.genes, m.size)

    m.genes = ['gene1','gene2','gene3','gene4']
    print(m.genes, m.size)

    print(m.cov)
