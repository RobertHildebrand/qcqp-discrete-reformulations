from glob import glob as glb
from os import path
from gurobipy import *
from Model_qcqp import Model_qcqp
from math import log, exp, isnan
from operator import mul
from functools import reduce
import json
import pprint
from Gurobi_2_Pyomo import *
import pyomo.environ as po

pp = pprint.PrettyPrinter(indent=2)

fsplit = path.splitext

### File settings; user input ###
fpath = sys.argv[1]
do_callback = int(sys.argv[2])

# fp: file directory
# fname: name of file, without extension
# fext: read-in extension
(fp, fnameext) = path.split(fpath)
(fname, fext) = fsplit(fnameext)

if len(fp) == 0:
    # If PWD is the folder with the instances, reference other folders as './../(stuff)' insted of '/../(stuff)'
    fp = '.'

# Make direct Gurobi filename
oname = re.sub('_L=(.*)', '', fname)  # Removes method-specific information
oname = re.sub('(.*)(_)(.*)', r'\1', oname)  # Removes last underscore and everything after it
opath = rf'{fp}/../original/{oname}.lp'

#print(oname)
#print(opath)

# Method and read/run/write settings
maxTime = 60*60*2  # Tell Lukas to add a 1 minute to the maximum runtime as a file IO allowance

rpath = rf'{fp}/../results/{fname}_{maxTime}s_{("no","with")[do_callback]}cb'

if len(fp) == 0:
    # If PWD is the folder with the instances, delete the slash at the start to prevent directories from root!
    rpath = rpath[1:]
    opath = opath[1:]

### No input from users past this point ###

# pname = 'minlp_ac_opf_nesta_case14_ieee__api_gurobi'
# fpath = rf'{fp}/{pname}{fext}'

# Read in necessary file(s)
mdl = read(fpath)
do_minimize = mdl.ModelSense == GRB.MINIMIZE

feas_sols = []
feas_objs = []
if do_callback:
    m = read(opath)
    # Set all variables to continuous, to circumvent something weird with Gurobi reading in variables as integers
    #   if there is no 'Subject To' line in the .lp file.
    vs = m.getVars()
    #for i in range(len(vs)):
    #    m.setAttr("vtype", [vs[i]], [GRB.CONTINUOUS])

    m.update()

    # Make pyomo version for use in callback
    G2P = Gurobi_2_Pyomo(m)
    G2P.convert_model()
    G2P.freeze_integers()

    m_po = G2P.mdl
    name2x = G2P.name2x
    onames = list(name2x.keys())
    obj = G2P.obj

    best_sol = 1e20 * (1 if do_minimize else -1)

    TC = po.TerminationCondition
    optmlTcons = [TC.globallyOptimal, TC.locallyOptimal, TC.optimal]


    def my_callback(model: Model, where):
        global best_sol
        if where == GRB.Callback.MIPSOL:
            xsol = model.cbGetSolution(orig_vars)
            for (ii, name) in enumerate(onames):
                name2x[name].value = xsol[ii]

            opt = po.SolverFactory('ipopt')
            results = opt.solve(m_po)

            if results.solver.status == po.SolverStatus.ok:

                oval = po.value(obj)
                sol_improved = (do_minimize and oval < best_sol or
                                not do_minimize and oval > best_sol)
                if sol_improved:
                    best_sol = oval
                    feas_sols.append({name: name2x[name].value for name in onames})
                    feas_objs.append(oval)
                    print(f'Full-problem solution found, obj={oval}')

if do_callback:
    orig_vars = [mdl.getVarByName(name) for name in name2x]

#if mthd == DIRECT:
#    mdl.setParam('NonConvex', 2)

mdl.setParam("LogFile", f'{rpath}.log')
mdl.setParam('TimeLimit', maxTime)
mdl.setParam('OutputFlag', 0)

if do_callback:
    mdl.setParam('MIPfocus', 1)
    mdl.optimize(callback=my_callback)
    #pp.pprint(feas_objs)
else:
    mdl.optimize()



if mdl.SolCount >= 1:
    mdl.write(f'{rpath}.mst')
    mdl.write(f'{rpath}.sol')
    bestSol = mdl.getObjective().getValue()
else:
    bestSol = float('inf')*(1 if do_minimize else -1)

if do_callback and len(feas_objs) >= 1:
    bestQCQPSol = (min(feas_objs) if do_minimize else max(feas_objs))
else:
    bestQCQPSol = float('inf')*(1 if do_minimize else -1)

rtime = mdl.getAttr(GRB.Attr.Runtime)
bestBnd = mdl.getAttr(GRB.Attr.ObjBound)
nodeCount = mdl.getAttr(GRB.Attr.NodeCount)

expDict = {'time': rtime, 'QCQPsol': bestQCQPSol, 'MIPsol': bestSol, 'bnd': bestBnd, 'nodeCount': nodeCount,
           'allQCQPsols': feas_objs}

with open(f'{rpath}.result', 'wt') as out:
    pprint.pprint(expDict, stream=out)

if do_callback:
    with open(f'{rpath}_fsols.cbsol', 'wt') as out:
        pprint.pprint(feas_sols, stream=out)

