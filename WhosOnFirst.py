from gurobipy import *

m = Model()

x = m.addVar(name='x')

m.addConstr(x <= 0)

m.setParam('OutputFlag', 0)
m.setParam('LogFile', 'idiocy.log')

m.optimize()
