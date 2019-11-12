# This is an abstract class for startegy

class BasicStartegy ():

    def __init__(self, condition, score_function,components_folder,Process_seed, contraction_allowed):
        self.condition = condition
        self.score_function = score_function
        self.components_folder =components_folder
        self.contraction_allowed = contraction_allowed
        self.Process_seed = Process_seed


    # an abstract class, that can be overriden
    # the details of the growing a module should be implemented in this function
    def __call__(self, current_module,externel_nodes, usednodes,prev_modules, fnc=None):
        pass
