from .basic import BasicStartegy
from  config import *

class OverlapStrategy(BasicStartegy):
    def __init__(self,condition, score_function,components_folder,
                    Process_seed = False,
                    contraction_allowed =True,
                    d_threshold=1):

        super(OverlapStrategy, self).__init__(condition, score_function,components_folder,Process_seed, contraction_allowed)
        self.d_threshold = d_threshold

    # Override
    def __call__(self,current_module,externel_nodes,articulation_points,
                    usednodes, degrees, G, used_as_seed ,
                    prev_modules,fnc= None):
        best_genes = []
        best_score = base_score = self.score_function(current_module)
        # which action to perform
        should_add = True
        # flag for stopping the process if there is no improvement in the score
        should_stop = False
        self.condition.value = current_module

        # Consider adding a gene
        for node in externel_nodes:

            n_used_times = 0
            module = current_module.copy().add_gene(node)
            try:
                # if the node not in usednodes, an error will be raised, meaning the node was not used
                # and except will hundle the error, otherwise continue as the node is already used
                node_module = usednodes[node]
                n_used_times = len(usednodes[node])

                mean_deg = np.array(degrees[node]).mean()
                deg_ = len(set(list(G.adj[node])).intersection(set(module.genes)))

                if deg_*self.d_threshold < mean_deg:
                    continue

            except:
                pass

            if not self.condition(module,fnc): continue
            score = self.score_function(module)

            if score > best_score:
                best_score = score
                best_genes=[node]
            elif score == best_score:
                best_genes.append(node)
            else: pass

        if self.contraction_allowed and current_module.size >1:

            # Consider removing a gene from the module
            for node in current_module.genes:
                # print('consider removing node: ', node)
                if (self.Process_seed and current_module.is_seed(node)): continue
                module = current_module.copy().remove_gene(node)
                score = self.score_function(module)
                if (score > base_score + 1e-12) or (score < base_score + 1e-12): continue

                # The following condition is necessary to avoid cases when a
				# tree-like module becomes disconnected due to the removal
				# of a non-leaf gene

                if node in articulation_points: continue

                if score > best_score:
                    should_add = False
                    best_score = score
                    best_genes=[node]
                elif score == best_score:
                    best_genes.append(node)
                else: pass


        if len(best_genes) == 0:
            # There was no gene that could imporve the score, so the process should stop
            should_stop = True

        return best_genes, should_add, should_stop




    def post_process_modules(self, modules, MAX_GENES, results_range):
        if results_range !=None:
            range_ = results_range

        elif MAX_GENES%100 != 0:
            range_ = list(range(100,MAX_GENES, 100))+[MAX_GENES]
        else:
            range_ = list(range(100,MAX_GENES+100, 100))

        for n in range_:
            tmp = []
            i=0
            while(i < n):
                for m in modules:
                    if n-i > m.size:
                        tmp.append(m.genes)
                        i += m.size
                    else:
                        tmp.append(m.genes[:n-i])
                        i += n-i
                        break
            
            with open(self.components_folder+'cc_n{}_k{}.txt'.format(n,k), 'w') as f:
                f.write('\n'.join([' '.join(m) for m in tmp]))
