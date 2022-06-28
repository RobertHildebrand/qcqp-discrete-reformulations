from gurobipy import *

fpath = r'D:\Users\Ben Beach\Google Drive\Research_grad\PHD Research\NN_QCQP\Lucas_collab\boxQP\LP-data\basic-LP\spar020-100-1_v2.lp'

mdl = read(fpath)

mdl.write(r'D:\Users\Ben Beach\Google Drive\Research_grad\PHD Research\NN_QCQP\Lucas_collab\boxQP\LP-data\basic-LP\spar020-100-1_v3.lp')

mdl.setParam('Nonconvex', 2)
mdl.optimize()


mdl = Model()
x = mdl.addVar(lb=0, ub=1, name='x')
y = mdl.addVar(lb=0, ub=1, name='y')
mdl.setObjective(-x**2 + y**2 - x*y + x - y)

mdl.write(r'D:\Users\Ben Beach\Google Drive\Research_grad\PHD Research\NN_QCQP\Lucas_collab\boxQP\LP-data\basic-LP\spar020-100-1_vtest.lp')