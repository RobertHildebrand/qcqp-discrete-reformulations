from glob import glob as glb
from os import path
from gurobipy import *
from Model_qcqp import Model_qcqp
from math import log, exp
import pprint

pp = pprint.PrettyPrinter(indent=2)

fsplit = path.splitext

# Flags from Model_qcqp
UNIVAR = Model_qcqp.UNIVAR
NMDT = Model_qcqp.NMDT
DNMDT = Model_qcqp.DNMDT
DIRECT = 9001

### File settings; user input ###
p_classes = ['basic', 'small', 'extended', 'extended2']
#p_classes = ['extended2']
out_path = r'./boxQP/mdl-data'


### IMPORTANT: name for run for result summary saving. User-specified so as not to become unweildy. ###
runstr = 'L3_ldas_0-1_20min'

### Model settings; user input ###
# All methods
L = 3

# DNMDT #
# ldas: list of lambda values to use for DNMT.
ldas = [0, 1]

# NMDT #
disc_left = 1
disc_right = 0

# Sawtooth #
L1 = 10
st_do_approx = 0

# List of methods #
mthds = [UNIVAR, NMDT, DNMDT, DIRECT]

mdict = {mthd: {} for mthd in mthds}

### Other setup ###
ldas_str = f'ldas={"-".join(str(i) for i in ldas)}'

dnmdt_str = f'_DNMDT_L={L}_ldas={ldas_str}'
nmdt_str = f'_NMDT_L={L}_discleft={disc_left}_discright={disc_right}'
univar_str = f'_UNIVAR_L={L}_L1={L1}_doapprox={st_do_approx}'
direct_str = f''

run_mdl = 1
maxTime = 1200


mthd_strs = {UNIVAR: univar_str, DNMDT: dnmdt_str, NMDT: nmdt_str, DIRECT: direct_str}
mthd_res_strs = {UNIVAR: univar_str[1:],
                 DNMDT: dnmdt_str[1:],
                 NMDT: nmdt_str[1:],
                 DIRECT: 'DIRECT'}

if not (disc_right or disc_left):
    disc_left = 0
for pcls in p_classes:
    # file path per pcls; user input
    myfpath = rf'./boxQP/LP-data/{pcls}-LP/*.lp'

    ### No input from users past this point ###
    fpaths = glb(myfpath)

    for fpath in fpaths:
        (fp, fnameext) = path.split(fpath)
        (pname, fext) = fsplit(fnameext)
        # fp: original file path
        # fname: name of file, without extension
        # fext: read-in extension

        m = read(fpath)

        # Set all variables to continuous, to circumvent something weird with Gurobi reading in variables as integers
        #   if there is no 'Subject To' line in the .lp file.
        vs = m.getVars()
        for i in range(len(vs)):
            m.setAttr("vtype", [vs[i]], [GRB.CONTINUOUS])

        m.update()

        mpath_exists = 1
        # Write all models for the current problem to files
        for mthd in mthds:
            mthd_str = mthd_strs[mthd]
            mpath = rf'{out_path}/{pcls}/{pname}{mthd_str}.lp'
            rpath = rf'{out_path}/{pcls}/Results/{pname}{mthd_str}'
            if mthd != DIRECT:
                if not(path.exists(mpath)):
                    mpath_exists = 0

                    myL = (2*L if mthd == NMDT else L)
                    mdl_holder = Model_qcqp(m, L=myL, L1=L1, ldas=[0], method=mthd,
                                            disc_right=disc_right, disc_left=disc_left, st_do_approx=st_do_approx)

                    mdl = mdl_holder.mdl
                    mdl.update()
                else:
                    mdl = read(mpath)

            else:
                mpath = rf'{out_path}/{pcls}/{pname}{mthd_str}.lp'
                if not(path.exists(mpath)):
                    mpath_exists = 0
                mdl = m


            if not mpath_exists:
                # write model to file
                mdl.write(rf'{out_path}/{pcls}/{pname}{mthd_str}.lp')

            if run_mdl:
                if mthd == DIRECT:
                    mdl.setParam('NonConvex', 2)
                mdl.setParam("LogFile", f'{rpath}.log')
                mdl.setParam('TimeLimit', maxTime)
                mdl.setParam('OutputFlag', 1)
                mdl.optimize()
                if mdl.SolCount >= 1:
                    mdl.write(f'{rpath}.mst')
                    mdl.write(f'{rpath}.sol')
                    bestSol = mdl.getObjective().getValue()
                else:
                    bestSol = (-1)**(mdl.ModelSense == GRB.MAXIMIZE)*float('inf')

                rtime = mdl.getAttr(GRB.Attr.Runtime)
                bestBnd = mdl.getAttr(GRB.Attr.ObjBound)
                nodeCount = mdl.getAttr(GRB.Attr.NodeCount)

                expDisc = (f'{pcls}_{pname}{mthd_str}:\t time={rtime},'
                           f'\t bestSol={bestSol},\t bestBound = {bestBnd}')

                expDict = {'time': rtime, 'sol': bestSol, 'bnd': bestBnd, 'nodeCount': nodeCount, 'expDisc': expDisc}

                with open(f'{rpath}_res.txt', 'wt') as out:
                    pprint.pprint(expDict, stream=out)

                mdict[mthd].update({f'{pcls}_{pname}': expDict})
            # During results processing, one will need to detect file nonexistence, then report a solution-free timeout.

times = {mthd: [] for mthd in mthds}
sols = {mthd: [] for mthd in mthds}
bnds = {mthd: [] for mthd in mthds}
gaps = {mthd: [] for mthd in mthds}
for mthd in mthds:
    rdict = mdict[mthd]

    times[mthd] = [rdict[key]['time'] for key in rdict]
    sols[mthd] = [rdict[key]['sol'] for key in rdict]
    bnds[mthd] = [rdict[key]['bnd'] for key in rdict]

    nexp = len(sols[mthd])
    I = range(nexp)

    # Reported relative gaps, not final relative gaps
    gaps[mthd] = [abs(sols[mthd][i]-bnds[mthd][i])/abs(sols[mthd][i])*100
                  if abs(sols[mthd][i])>=1e-10 else float('inf')
                  for i in I]

    rstr = rf'{out_path}/{mthd_res_strs[mthd]}'
    with open(rf'{rstr}_allres.txt', 'wt') as out:
        pprint.pprint(rdict, stream=out)


BB = [max(bnds[mthd][i] for mthd in mthds) for i in I]
TO = {(mthd, i): (times[mthd][i] >= maxTime) for mthd in mthds for i in I}

I_Solved = [i for i in I if all(not TO[mthd, i] for mthd in mthds)]
I_Contested = [i for i in I if any(TO[mthd, i] for mthd in mthds) and not all(TO[mthd, i] for mthd in mthds)]
I_Unsolved = [i for i in I if all(TO[mthd, i] for mthd in mthds)]

SLV = 0
CNT = 1
USLV = 2

clss = [SLV, CNT, USLV]

Icls = {SLV: I_Solved, CNT: I_Contested, USLV: I_Unsolved}
clsStr = {SLV: 'solved', CNT: 'contested', USLV: 'unsolved'}

for cls in clss:
    Is = Icls[cls]

    tol = 1e-10
    n_cls = len(Is)

    if n_cls >= 1:
        BB_cls = {mthd_res_strs[mthd]: sum(abs(bnds[mthd][i] - BB[i]) <= 1e-10 for i in Is)
                  for mthd in mthds}
        TO_cls = {mthd_res_strs[mthd]: sum(times[mthd][i] >= maxTime for i in Is)
                  for mthd in mthds}

        tshift = min(times[mthd][i] for i in Is for mthd in mthds)

        # Compute average time via shifted geometric mean
        tave_cls = {mthd_res_strs[mthd]: exp(sum(log(times[mthd][i]-tshift+1) for i in Is)/n_cls)+tshift-1
                    for mthd in mthds}
    else:
        BB_cls = {mthd_res_strs[mthd]: 0
                  for mthd in mthds}
        TO_cls = {mthd_res_strs[mthd]: 0
                  for mthd in mthds}
        tave_cls = {mthd_res_strs[mthd]: -1
                    for mthd in mthds}

    rdict = {'n': n_cls, 'BB': BB_cls, 'TO': TO_cls, 'time': tave_cls}

    rstr = rf'{out_path}/res_{runstr}_{clsStr[cls]}.txt'
    with open(rstr, 'wt') as out:
        pprint.pprint(rdict, stream=out)

