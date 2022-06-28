from gurobipy import *

# The formulations in this module allow the user to provide the auxiliary variables. This give the user the option
#   To use the same auxiliary variables among different terms.
class Model_xsq:
    def __init__(self, mdl, x, y, L, lb=0, ub=1, alpha=None, g=None):
        self.mdl = mdl
        self.x = x
        self.y = y
        self.L = L
        self.lb = lb
        self.ub = ub
        self.Lset = range(L)

        if alpha is not None:
            self.alpha = alpha
        if g is not None:
            self.g = g

    def sawtooth(self, mdl, x, y, lb=0, ub=1):
        # Apply sawtooth formulation to model y=x^2
        pass

    def sawtooth_ip(self, mdl, x, y, lb=0, ub=3):
        # Apply the IP-sawtooth formulation to model y=x^2, x in [lb,ub], x integer
        pass

    def sawtooth_ip_sym(self, mdl, x, y, lb=0, ub=3):
        # Apply the symmetric version of the IP-sawtooth formulation to model y=x^2, x in [lb,ub], x integer
        pass

    def nmdt(self, ):
        # Apply NMDT to model y=x^2
        pass

    def dnmdt(self, ):
        # Apply DNMDT to model y=x^2
        pass

    def set(self, L=None, lb=None, ub=None):
        if L is not None:
            self.L=L
        if lb is not None:
            self.lb=lb
        if ub is not None:
            self.ub=ub