from .basic import BasicStartegy
from  config import *

class OverlapStrategy(BasicStartegy):
    def __init__(self,condition, score_function,components_folder, Process_seed = False, contraction_allowed =True, d_threshold=1):
        super(OverlapStrategy, self).__init__(condition, score_function,components_folder,Process_seed, contraction_allowed)
        self.d_threshold = d_threshold

    # Override
    def __call__(self,current_module,externel_nodes,articulation_points, usednodes, degrees, G, used_as_seed ,prev_modules,fnc= None):
        best_genes = []
        best_score = base_score = self.score_function(current_module)
        # which action to perform
        should_add = True
        # flag for stopping the process if there is no improvement in the score
        should_stop = False
        self.condition.value = current_module

        # Consider adding a gene
        for node in externel_nodes:
            # try:
            #     # if the node was used as a seed, don't consider it
            #     _ = used_as_seed[node]
            #     continue
            # except:
            #     pass
            n_used_times = 0
            module = current_module.copy().add_gene(node)
            try:
                # if the node not in usednodes, an error will be raised, meaning the node was not used
                # and except will hundle the error, otherwise continue as the node is already used
                node_module = usednodes[node]
                n_used_times = len(usednodes[node])
                # f_overlap = False
                # If a node (gene) has been used N times before, don't consider it
                # N =10
                # if len(node_module) > N:
                #     continue
                # for idx in node_module:
                #     if not module.overlap(prev_modules[idx], overlap_threshold):
                #         # print('Adding gene {} increase the overlap with module {}'.format(node, idx))
                #         f_overlap = True
                #         break
                # This else belongs to the for loop
                #
                # else:

                mean_deg = np.array(degrees[node]).mean()
                deg_ = len(set(list(G.adj[node])).intersection(set(module.genes)))

                # if node == 'TP53':
                #     print('#used {}, deg_ {}, mean_deg {}'.format(n_used_times, deg_, mean_deg))
                if deg_*self.d_threshold < mean_deg:
                    continue


                if f_overlap:
                    # This continue for the outerloop line 19. it means check another gene
                    continue

            except:
                pass
            # t = 0.95 + n_used_times*5
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

        # print('Addition is considered') if should_add else  print('Removing is considered')
        if len(best_genes) == 0:
            # There was no gene that could imporve the score, so the process should stop
            should_stop = True

        return best_genes, should_add, should_stop

    # This method was copied from the Nonoverlapp strategy, where we only care
    # about the total genes (sum of modules sizes), regarless of the uniquness.
    #  It should be modified for proper running
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
            # if n == 100:
                # print(tmp)
            with open(self.components_folder+'cc_n{}_k{}.txt'.format(n,k), 'w') as f:
                f.write('\n'.join([' '.join(m) for m in tmp]))

    # def post_process_modules(self, modules, MAX_GENES):
    #     unique_genes = set()
    #     if MAX_GENES%100 != 0:
    #         range_ = list(range(100,MAX_GENES, 100))+[MAX_GENES]
    #     else:
    #         range_ = list(range(100,MAX_GENES+100, 100))
    #
    #     for n in range_:
    #         tmp = []
    #         i=0
    #         while(i < n):
    #             for m in modules:
    #                 # print('i: ', i)
    #                 tmp_unique = unique_genes.union(set(m.genes))
    #                 if len(tmp_unique) <= n:
    #                     tmp.append(m.genes)
    #                     unique_genes = unique_genes.union(set(m.genes))
    #                     i = len(unique_genes)
    #                 else:
    #                     tmp_unique = unique_genes
    #                     t = []
    #                     for g in m.genes:
    #                         if len(tmp_unique.union(set([g]))) <= n:
    #                             tmp_unique = tmp_unique.union(set([g]))
    #                             t.append(g)
    #                             i+=1
    #
    #                     tmp.append(t)
    #                     break
    #         if n == 100:
    #             print(tmp)
    #         with open(self.components_folder+'cc_n{}_k{}.txt'.format(n,k), 'w') as f:
    #             f.write('\n'.join([' '.join(m) for m in tmp]))
