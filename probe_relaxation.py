from glob import glob as glb
from os import path
from gurobipy import *
from Model_qcqp import Model_qcqp
import pprint

pp = pprint.PrettyPrinter(indent=2)

UNIVAR = Model_qcqp.UNIVAR
NMDT = Model_qcqp.NMDT
DNMDT = Model_qcqp.DNMDT
BIN2 = Model_qcqp.BIN2
BIN3 = Model_qcqp.BIN3
MCCORMICK = Model_qcqp.MCCORMICK
DIRECT = 9001


L = 0
L1 = 0
xl = -3
xu = 5
yl = -3
yu = 5

lx = xu-xl
ly = yu-yl


m = Model()

x = m.addVar(lb=xl, ub=xu, name='x')
y = m.addVar(lb=yl, ub=yu, name='y')
z = m.addVar(lb=-GRB.INFINITY, name='z')

m.setObjective(z)
m.addConstr(z == x*y, name='mult')
#m.addConstr(x*(x-y) >= 2, name='m')
m.setParam('NonConvex', 2)
m.optimize()

#methods = [BIN2, BIN3, UNIVAR, NMDT, DNMDT, MCCORMICK]
methods = [UNIVAR, DNMDT]


mdls = {method: Model_qcqp(m, L=L, L1=L1, ldas=[0, 1], method=method, disc_right=0, st_do_approx=0).mdl
        for method in methods}

nx = 5
ny = 5
xldas = [i/(nx-1) for i in range(nx)]
yldas = [i/(ny-1) for i in range(ny)]
xs = [lda*xl + (1-lda)*xu for lda in xldas]
ys = [lda*yl + (1-lda)*yu for lda in yldas]
#xs = [.4, .8]
#ys = [.4, .8]

lbs = {method: {(a, b): 0 for a in xs for b in ys}
       for method in methods}

ubs = {method: {(a, b): 0 for a in xs for b in ys}
       for method in methods}

for method in methods:
    mdl = mdls[method]
    mdl.setParam('OutputFlag', 0)
    mdl.update()

    mvs = mdl.getVars()
    print(mvs)
    x = [v for v in mvs if v.varName == 'x'][0]
    y = [v for v in mvs if v.varName == 'y'][0]
    z = [v for v in mvs if v.varName == 'z'][0]
    for a in xs:
        for b in ys:
            if a==b and a==2:
                mdl.setParam('OutputFlag', 1)
            else:
                mdl.setParam('OutputFlag', 0)

            tc1 = mdl.addConstr(x == a)
            tc2 = mdl.addConstr(y == b)

            mdl.setObjective(z)
            mdl.optimize()
            #lbs[method][a,b] = mdl.getObjective().getValue()
            if mdl.Status == GRB.OPTIMAL:
                lbs[method][a, b] = mdl.getObjective().getValue()
            else:
                lbs[method][a, b] = -9001


            mdl.setObjective(-z)
            mdl.optimize()
            if mdl.Status == GRB.OPTIMAL:
                ubs[method][a, b] = -mdl.getObjective().getValue()
            else:
                ubs[method][a, b] = 9001

            mdl.remove(tc1)
            mdl.remove(tc2)

lb_mstrs = {(a, b): [f'{round(lbs[method][a,b],6)}' for method in methods]
            for a in xs for b in ys}
ub_mstrs = {(a, b): [f'{round(ubs[method][a,b],6)}' for method in methods]
            for a in xs for b in ys}
pp.pprint({(a, b): ', '.join(lb_mstrs[a, b])
           for a in xs for b in ys if len(set(lb_mstrs[a, b])) >= 2})
pp.pprint({(a, b): ', '.join(ub_mstrs[a, b])
           for a in xs for b in ys if len(set(ub_mstrs[a, b])) >= 2})