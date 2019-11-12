# This is the class that contains the condition for considering a node (gene) for adding or removing
# it has a function condition returns True or False
# it has a set_condition function
from .functions import *
import numpy as np

class Condition():

    def __init__(self, condition_type='cov/mutex',threshold = 0.9):

        if not isinstance(condition_type, str):
            raise Exception('The condition_type should be of type str')
        if not condition_type in fncs.keys():
            raise Exception('The condition_type should be one of the following: {}'.format(list(functions.keys())))
        self.function = fncs[condition_type]
        self.threshold = threshold
        self.condition_type = condition_type
        self._condition = None

    @property
    def value(self):
        return self._condition

    @value.setter
    def value(self,current_module):
        self._condition  = self.function(current_module)


    def __call__(self, considered_module,threshold=None, fnc =None):
        if not threshold == None:
            self.threshold = threshold

        if fnc:
            return fnc(considered_module, self.threshold) > self._condition
        else:
            return self.function(considered_module, self.threshold) > self._condition
