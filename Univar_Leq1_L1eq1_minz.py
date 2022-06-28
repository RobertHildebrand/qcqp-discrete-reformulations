from gurobipy import *

m = Model()

ly = 0

y = m.addVar(lb=ly, ub=ly+1/2)

zy = m.addVar(lb=-GRB.INFINITY)
zp1 = m.addVar(lb=0)
z = m.addVar(lb=-GRB.INFINITY, obj=1)

if ly == 0:
    m.addConstr(zy <= y/2)
else:
    m.addConstr(zy <= (3*y-1)/2)

m.addConstr(zp1 >= 3*y-2-1/4)
m.addConstr(zp1 >= 2*y-1)
m.addConstr(zp1 >= y-1/4)
#m.addConstr(zp1 >= 4*y-4)
m.addConstr(z >= zp1 - zy)

m.optimize()

print(f'y={y.x},\nzy={zy.x},\nzp1={zp1.x},\nz={z.x}')

