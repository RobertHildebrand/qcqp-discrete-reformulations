from glob import glob as glb
from os import path
from gurobipy import *
from Model_qcqp import Model_qcqp
from math import log, exp, isnan
from operator import mul
from functools import reduce
import json
import pprint

pp = pprint.PrettyPrinter(indent=2)

fsplit = path.splitext

# Flags from Model_qcqp
UNIVAR = Model_qcqp.UNIVAR
NMDT = Model_qcqp.NMDT
DNMDT = Model_qcqp.DNMDT
BIN2 = Model_qcqp.BIN2
BIN3 = Model_qcqp.BIN3
MCCORMICK = Model_qcqp.MCCORMICK
DIRECT = 9001

### File settings; user input ###
out_path = r'./ACOPF/mdl-data'


### IMPORTANT: name for run for result summary saving. User-specified so as not to become unweildy. ###

### Model settings; user input ###
# All methods
L = 1

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
run_mthds = []
read_mthds = [BIN2, BIN3, UNIVAR, NMDT, DNMDT]
#run_mthds = [BIN2, BIN3, UNIVAR, NMDT, DNMDT]
#read_mthds = []

#mthds = [BIN2, BIN3, UNIVAR, NMDT, DNMDT, MCCORMICK]
mthds = run_mthds + read_mthds
#mthds = [MCCORMICK, NMDT, DNMDT]
#mthds = [BIN2, BIN3, UNIVAR]
#mthds = [MCCORMICK, BIN2]

mdict = {mthd: {} for mthd in mthds}

def long2short(longstr):
    return re.sub('^minlp_ac_opf_|_gurobi(.*)', '', longstr)

with open(f'{out_path}/best_sols.txt', "r") as file:
    best_sols = eval(file.read())


### Other setup ###
ldas_str = f'ldas={"-".join(str(i) for i in ldas)}'

dnmdt_str = f'_DNMDT_L={L}_ldas={ldas_str}'
nmdt_str = f'_NMDT_L={L}_discleft={disc_left}_discright={disc_right}'
univar_str = f'_UNIVAR_L={L}_L1={L1}_doapprox={st_do_approx}'
bin2_str = f'_BIN2_L={L}_L1={L1}'
bin3_str = f'_BIN3_L={L}_L1={L1}'
mc_str = f'_MCCORMICK'
direct_str = f''

run_mdl = 1
#maxTime = 1200
maxTime = 600

runstr = f'L{L}_ldas_0-1_{maxTime}s'

mthd_strs = {UNIVAR: univar_str, DNMDT: dnmdt_str, NMDT: nmdt_str, DIRECT: direct_str,
             BIN2: bin2_str, BIN3: bin3_str, MCCORMICK: mc_str}
mthd_res_strs = {UNIVAR: univar_str[1:],
                 DNMDT: dnmdt_str[1:],
                 NMDT: nmdt_str[1:],
                 BIN2: bin2_str[1:],
                 BIN3: bin3_str[1:],
                 MCCORMICK: 'MCCORMICK',
                 DIRECT: 'DIRECT'}

if not (disc_right or disc_left):
    disc_left = 0

# File path; user input
myfpath = r'./ACOPF/lpfiles/*.lp'

### No input from users past this point ###
fpaths = glb(myfpath)

for fpath in fpaths:
    (fp, fnameext) = path.split(fpath)
    (pname, fext) = fsplit(fnameext)
    # fp: original file path
    # fname: name of file, without extension
    # fext: read-in extension

    if len(run_mthds) >= 1:
        m = read(fpath)
        # Set all variables to continuous, to circumvent something weird with Gurobi reading in variables as integers
        #   if there is no 'Subject To' line in the .lp file.
        vs = m.getVars()
        #for i in range(len(vs)):
        #    m.setAttr("vtype", [vs[i]], [GRB.CONTINUOUS])

        m.update()

    remodel_anyway = 0  # Bypass to re-create all models, in case of reruns with caught bugs.
    mpath_exists = 1
    # Write all models for the current problem to files
    for mthd in mthds:
        mthd_str = mthd_strs[mthd]
        mpath = rf'{out_path}/{pname}{mthd_str}.lp'
        rpath = rf'{out_path}/Results/{pname}{mthd_str}'
        if mthd in run_mthds:
            if mthd != DIRECT:
                if not(path.exists(mpath)) or remodel_anyway:
                    mpath_exists = 0

                    myL = (2*L if mthd == NMDT else L)
                    mdl_holder = Model_qcqp(m, L=myL, L1=L1, ldas=[0], method=mthd,
                                            disc_right=disc_right, disc_left=disc_left, st_do_approx=st_do_approx)

                    mdl = mdl_holder.mdl
                    mdl.update()
                else:
                    mdl = read(mpath)

            else:
                mpath = rf'{out_path}/{pname}{mthd_str}.lp'
                if not(path.exists(mpath)):
                    mpath_exists = 0
                mdl = m

            # 1-space looked weird here, but the Pycharm linter was complaining, so I added this comment.
            if not mpath_exists or remodel_anyway:
                # write model to file
                mdl.write(rf'{out_path}/{pname}{mthd_str}.lp')

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

                expDisc = (f'{pname}{mthd_str}:\t time={rtime},'
                           f'\t bestSol={bestSol},\t bestBound = {bestBnd}')

                expDict = {'time': rtime, 'sol': bestSol, 'bnd': bestBnd, 'nodeCount': nodeCount, 'expDisc': expDisc}

                with open(f'{rpath}_{maxTime}s_res.txt', 'wt') as out:
                    pprint.pprint(expDict, stream=out)

                mdict[mthd].update({f'{pname}': expDict})
        else:
            inf = float('inf')
            with open(f'{rpath}_{maxTime}s_res.txt', "r") as file:
                expdict = eval(file.read())
            mdict[mthd].update({f'{pname}': expdict})
        # During results processing, one will need to detect file nonexistence, then report a solution-free timeout.

times = {mthd: [] for mthd in mthds}
sols = {mthd: [] for mthd in mthds}
bnds = {mthd: [] for mthd in mthds}
gaps = {mthd: [] for mthd in mthds}

I = []
for mthd in mthds:
    rdict = mdict[mthd]

    bsols = [float(best_sols[long2short(key)]) for key in rdict]

    times[mthd] = [rdict[key]['time'] for key in rdict]
    sols[mthd] = [rdict[key]['sol'] for key in rdict]
    bnds[mthd] = [rdict[key]['bnd'] for key in rdict]

    nexp = len(sols[mthd])
    I = range(nexp)

    # Reported relative gaps, wrt best known solution to the problem, as provided by Robert Burlacu
    gaps[mthd] = [abs(bsols[i]-bnds[mthd][i])/abs(bsols[i])*100
                  if abs(bsols[i]) >= 1e-10 else float('inf')
                  for i in I]

    rstr = rf'{out_path}/{mthd_res_strs[mthd]}'
    with open(rf'{rstr}_{maxTime}s_allres.txt', 'wt') as out:
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
        nsols_cls = {mthd_res_strs[mthd]: len([sols[mthd][i] for i in Is if sols[mthd][i] < float('inf')])
                     for mthd in mthds}

        shift = max(1e-4, min(gaps[mthd][i] for i in Is for mthd in mthds))

        Is_gap = [i for i in Is if not any(isnan(gaps[mthd][i]) for mthd in mthds)]
        #gap_sgms = {mthd_res_strs[mthd]: reduce(mul, [gaps[mthd][i]+shift for i in Is], 1)**(1/len(Is))-shift
        #            for mthd in mthds}
        gap_sgms = {mthd_res_strs[mthd]: exp(sum(log(gaps[mthd][i]+shift) for i in Is_gap)/len(Is_gap))-shift
                    for mthd in mthds}

        bsols = [float(best_sols[long2short(key)]) for key in mdict[mthds[0]]]
    else:
        BB_cls = {mthd_res_strs[mthd]: 0
                  for mthd in mthds}
        TO_cls = {mthd_res_strs[mthd]: 0
                  for mthd in mthds}
        tave_cls = {mthd_res_strs[mthd]: -1
                    for mthd in mthds}
        nsols_cls = {mthd_res_strs[mthd]: -1
                     for mthd in mthds}
        gap_sgms = {mthd_res_strs[mthd]: -1
                    for mthd in mthds}

    rdict = {'n': n_cls, 'BB': BB_cls, 'TO': TO_cls, 'time': tave_cls, 'nsols': nsols_cls, 'gap_sgms': gap_sgms}
    rstr = rf'{out_path}/res_{runstr}_{maxTime}s_{clsStr[cls]}.txt'
    with open(rstr, 'wt') as out:
        pprint.pprint(rdict, stream=out)


