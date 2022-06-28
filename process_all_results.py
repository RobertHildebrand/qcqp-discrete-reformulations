from glob import glob as glb
from os import path
from gurobipy import *
from Model_qcqp import Model_qcqp
from math import log, exp, isnan
from operator import mul
from functools import reduce
import json
import pprint
from os.path import exists
import pandas
from numpy import nanmax, nanmin, maximum

pp = pprint.PrettyPrinter(indent=2)

fsplit = path.splitext

# Flags from Model_qcqp
UNIVAR = Model_qcqp.UNIVAR
NMDT = Model_qcqp.NMDT
DNMDT = Model_qcqp.DNMDT
TNMDT = Model_qcqp.TNMDT
TDNMDT = Model_qcqp.TDNMDT
BIN2 = Model_qcqp.BIN2
BIN3 = Model_qcqp.BIN3
MCCORMICK = Model_qcqp.MCCORMICK
DIRECT = 9001

### File settings; user input ###
out_path = r'./reuslts_240322'


### IMPORTANT: name for run for result summary saving. User-specified so as not to become unweildy. ###

### Model settings; user input ###
# All methods
L = 6

# DNMDT #
# ldas: list of lambda values to use for DNMT.
ldas = [0, 1]

# NMDT #
disc_left = 1
disc_right = 0

# Sawtooth #
L1 = 10
st_do_approx = 0

maxTime = 14400
#maxTime = 7200
callback_included = True

acopf_only = False
boxqp_only = True

if acopf_only:
    boxqp_only = False

datastr = ('ACOPF_' if acopf_only else
           'boxQP_' if boxqp_only else
           '')

startstr = ('minlp_ac_opf_nesta' if acopf_only else
            'spar' if boxqp_only else
            '')

callback_str = ('nocb', 'withcb')[callback_included]

# List of methods #
mthds = [BIN2, BIN3, UNIVAR, NMDT, DNMDT, TNMDT, TDNMDT]

mdict = {mthd: {} for mthd in mthds}


# ACOPF-specific stuff
def long2short(longstr):
    return re.sub('^minlp_ac_opf_|_gurobi(.*)', '', longstr)


def long2short_boxqp(longstr):
    return re.sub('^(spar[^_]*)_.*', r'\1', longstr)


def best_bqp_istr(dotin_str):
    return re.sub('^([^-]*)-|.in$', '', dotin_str)


if acopf_only:
    with open(f'{out_path}/best_sols.txt', "r") as file:
        best_sols = eval(file.read())


if boxqp_only:
    boxqp_data = pandas.read_csv(f'{out_path}/boxQP_sawtooth_results.csv')

    agg_df = boxqp_data.groupby('instance').agg({'true_objective': nanmin, 'dual_bound': nanmax})
    agg_dict = agg_df.T.to_dict()

    bbnds_dict = {best_bqp_istr(key): agg_dict[key]['dual_bound']
                  for key in agg_dict}

    bsols_dict = {best_bqp_istr(key): agg_dict[key]['true_objective']
                  for key in agg_dict}


### Other setup ###
ldas_str = f'ldas={"-".join(str(i) for i in ldas)}'

mthd_name_strs = {
    DNMDT: 'DNMDT',
    NMDT: 'NMDT',
    TDNMDT: 'T-DNMDT',
    TNMDT: 'T-NMDT',
    UNIVAR: 'Univar',
    BIN2: 'Bin2',
    BIN3: 'Bin3',
    MCCORMICK: 'McC',
    DIRECT: 'GRB'
}

# NMTT: NMDT, with pure quadratic terms tightened with sawtooth LB improvement
mthd_strs = {
    DNMDT: f'_DNMDT_L={L}_ldas={ldas_str}',
    NMDT: f'_NMDT_L={L}_discleft={disc_left}_discright={disc_right}',
    TDNMDT: f'_T-DNMDT_L={L}_ldas={ldas_str}',
    TNMDT: f'_T-NMDT_L={L}_discleft={disc_left}_discright={disc_right}',
    UNIVAR: f'_UNIVAR_L={L}_L1={L1}_doapprox={st_do_approx}',
    BIN2: f'_BIN2_L={L}_L1={L1}',
    BIN3: f'_BIN3_L={L}_L1={L1}',
    MCCORMICK: f'_MCCORMICK',
    DIRECT: f''
}

mthd_res_strs = {mthd: (mthd_strs[mthd] if mthd in [MCCORMICK, DIRECT]
                        else mthd_strs[mthd][1:])
                 for mthd in mthds}

fstrs = {mthd: f'{mthd_strs[mthd]}_{maxTime}s_{callback_str}'
         for mthd in mthds}

#maxTime = 1200

runstr = f'{datastr}L{L}_L1{L1}_{maxTime}s_{callback_str}'


if not (disc_right or disc_left):
    disc_left = 0

fpath_str = {mthd: rf'./reuslts_240322/results/{startstr}*{fstrs[mthd]}.log'
             for mthd in mthds}

'''
fpath_str = {mthd: rf'./results_120122/results/*_{mthd_name[mthd]}*nocb.result'
             for mthd in mthds}
'''

mdict = {mthd: {} for mthd in mthds}

inf = float('inf')
for mthd in mthds:
    # Read data from files
    # File path; user input
    fpaths = glb(fpath_str[mthd])

    for fpath in fpaths:
        (fp, fnameext) = path.split(fpath)
        (pname, fext) = fsplit(fnameext)
        fpath = rf'{fp}/{pname}.result'

        if exists(fpath):
            with open(fpath, "r") as file:
                expdict = eval(file.read())
        else:
            print(f'Errored file: {fpath}')
            expdict = {
                'MIPsol': inf,
                'QCQPsol': inf,
                'allQCQPsols': [],
                'bnd': -inf,
                'nodeCount': 0,
                'time': maxTime
            }

        mdict[mthd][f'{pname}'] = expdict


times = {mthd: [] for mthd in mthds}
sols = {mthd: [] for mthd in mthds}
Qsols = {mthd: [] for mthd in mthds}
bnds = {mthd: [] for mthd in mthds}
gaps = {mthd: [] for mthd in mthds}

I = []

for mthd in mthds:
    rdict = mdict[mthd]

    if mthd == mthds[0]:
        if boxqp_only:
            inames = [long2short_boxqp(key) for key in rdict]
            bsols = [float(bsols_dict[long2short_boxqp(key)]) for key in rdict]
            bbnds = [float(bbnds_dict[long2short_boxqp(key)]) for key in rdict]

        if acopf_only:
            inames = [long2short(key) for key in rdict]
            bsols = [float(best_sols[long2short(key)]) for key in rdict]

    times[mthd] = [rdict[key]['time'] for key in rdict]
    sols[mthd] = [rdict[key]['MIPsol'] for key in rdict]
    bnds[mthd] = [rdict[key]['bnd'] for key in rdict]

    if mthd == mthds[0] and boxqp_only:
        imprs = {inames[i]: [bsols[i]] for (i, bnd) in enumerate(bnds[mthd])}

    if callback_included:
        Qsols[mthd] = [rdict[key]['QCQPsol'] for key in rdict]

    if boxqp_only:
        for (i, bnd) in enumerate(bnds[mthd]):
            pbnd = bbnds[i]
            if bnd not in (-inf, inf) and bnd > pbnd:
                if abs(bnd-pbnd) >= 100*abs(pbnd)*1e-8:
                    print(f'New best bound found for instance {inames[i]} using {mthd_name_strs[mthd]}.')
                    print(f'Improvement: {pbnd}->{bnd}, {(bnd-pbnd)/pbnd*100}%.')
                if callback_included:
                    imprs[inames[i]].append((mthd_name_strs[mthd], Qsols[mthd][i], (pbnd, bnd), (bnd-pbnd)/pbnd*100))
                else:
                    imprs[inames[i]].append((mthd_name_strs[mthd], bsols[i], (pbnd, bnd), (bnd-pbnd)/pbnd*100))
                bbnds[i] = bnd

    if callback_included:
        if acopf_only or boxqp_only:
            for (i, sol) in enumerate(Qsols[mthd]):
                psol = bsols[i]
                if sol not in (-inf, inf) and sol < bsols[i]:
                    if (psol-sol) >= abs(psol)*100*1e-8:
                        print(f'New best solution found for instance {inames[i]} using {mthd_name_strs[mthd]}.')
                        print(f'Improvement: {psol}->{sol}, {(psol-sol)/psol*100}%.')
                    bsols[i] = sol

    for i in range(len(rdict)):
        if bnds[mthd][i] == inf:
            bnds[mthd][i] = -inf

    nexp = len(sols[mthd])
    I = range(nexp)

    # Reported relative gaps, wrt best known solution to the problem, as provided by Robert Burlacu

    with open(rf'{out_path}/allres/{datastr}{fstrs[mthd]}_allres.txt', 'wt') as out:
        pprint.pprint(rdict, stream=out)

if boxqp_only:
    pp.pprint({key: imprs[key] for key in imprs if len(imprs[key]) >= 2})

for mthd in mthds:
    if acopf_only:
        gaps[mthd] = [abs(bsols[i]-bnds[mthd][i])/abs(bsols[i])*100
                      if abs(bsols[i]) >= 1e-10 else inf
                      for i in I]


if boxqp_only:
    #bsols = [min(Qsols[mthd][i] for mthd in mthds) for i in I]
    #bbnds = [max(bnds[mthd][i] for mthd in mthds) for i in I]
    if callback_included:
        sgaps = {}

    for mthd in mthds:
        gaps[mthd] = [abs(bsols[i]-bnds[mthd][i])/abs(bsols[i])*100
                      if abs(bsols[i]) >= 1e-10 else inf
                      for i in I]
        if callback_included:
            sgaps[mthd] = [abs(Qsols[mthd][i]-bbnds[i])/abs(Qsols[mthd][i])*100
                           if abs(Qsols[mthd][i]) >= 1e-10 else inf
                           for i in I]

if acopf_only and False:
    for i in I:
        if all(gaps[mthd][i] == inf for mthd in mthds):
            continue

        max_gap = max(gaps[mthd][i] for mthd in mthds if gaps[mthd][i] != inf)
        for mthd in mthds:
            if gaps[mthd][i] == inf:
                # Upon erroring, gap=2*max_gap
                gaps[mthd][i] = 2*max_gap

BB = [max(bnds[mthd][i] for mthd in mthds) for i in I]

rel_bnds = {mthd: [abs(bnds[mthd][i] - BB[i])/abs(BB[i]) for i in I]
            for mthd in mthds}

#pp.pprint([f'BB:{BB[i]}, {", ".join(f"{mthd_strs[mthd]}:{rel_bnds[mthd][i]}" for mthd in mthds)}' for i in I])


TO = {(mthd, i): (times[mthd][i] >= maxTime) for mthd in mthds for i in I}

I_Solved = [i for i in I if all(not TO[mthd, i] for mthd in mthds)]
I_Contested = [i for i in I if any(TO[mthd, i] for mthd in mthds) and not all(TO[mthd, i] for mthd in mthds)]
I_Unsolved = [i for i in I if all(TO[mthd, i] for mthd in mthds)]

SLV = 0
CNT = 1
USLV = 2
ALL = 3

#clss = [SLV, CNT, USLV, ALL]
clss = [SLV, CNT, USLV]

Icls = {SLV: I_Solved, CNT: I_Contested, USLV: I_Unsolved, ALL: I}
cls_str = {SLV: 'solved', CNT: 'contested', USLV: 'unsolved'}

tabstr = ''

print(Icls)

# ACOPF, L1, noCB classes
#Icls = {SLV:  [3, 5, 22, 23, 24, 25, 26, 27, 29, 33, 34, 35, 36, 37, 38, 42, 43, 44, 45, 46, 47, 48, 49, 56, 57, 58],
#        CNT:  [4, 9, 11, 12, 13, 14, 16, 21, 28, 39, 40, 41, 53, 54, 55],
#        USLV: [0, 1, 2, 6, 7, 8, 10, 15, 17, 18, 19, 20, 30, 31, 32, 50, 51, 52]}

# ACOPF, L2, noCB classes
#Icls = {SLV: [33, 34, 35, 36, 37, 38, 42, 43, 44, 45, 46, 47, 48, 49, 56, 57, 58],
#        CNT: [3, 4, 5, 21, 22, 23, 24, 25, 26, 27, 28, 29, 54],
#        USLV: [0, 1, 2, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 30, 31, 32, 39, 40, 41, 50, 51, 52, 53, 55]}


#pp.pprint([[rel_bnds[UNIVAR][i]] for i in I_Contested])

for cls in clss:
    Is = Icls[cls]

    tol = 1e-10
    n_cls = len(Is)

    if n_cls >= 1:
        BB_cls = {mthd_res_strs[mthd]: sum(rel_bnds[mthd][i] <= 1e-5 for i in Is)
                  for mthd in mthds}
        TO_cls = {mthd_res_strs[mthd]: sum(times[mthd][i] >= maxTime for i in Is)
                  for mthd in mthds}

        tshift = min(times[mthd][i] for i in Is for mthd in mthds)

        # Compute average time via shifted geometric mean
        tave_cls = {mthd_res_strs[mthd]: exp(sum(log(times[mthd][i]-tshift+1) for i in Is)/n_cls)+tshift-1
                    for mthd in mthds}
        nsols_cls = {mthd_res_strs[mthd]: len([sols[mthd][i] for i in Is if sols[mthd][i] < float('inf')])
                     for mthd in mthds}

        if callback_included:
            nQsols_cls = {mthd_res_strs[mthd]: len([Qsols[mthd][i] for i in Is if Qsols[mthd][i] < float('inf')])
                          for mthd in mthds}

        bgap_shift = 1e-5
        Is_gap = [i for i in Is if not any(isnan(rel_bnds[mthd][i]) or math.isinf(rel_bnds[mthd][i]) for mthd in mthds)]
        BB_gaps = {mthd_res_strs[mthd]: exp(sum(log(rel_bnds[mthd][i]+bgap_shift) for i in Is_gap)/len(Is_gap))-bgap_shift
                   for mthd in mthds}

        if acopf_only or boxqp_only:
            shift = max(1e-4, min(gaps[mthd][i] for i in Is for mthd in mthds))

            Is_gap = [i for i in Is if not any(isnan(gaps[mthd][i]) or gaps[mthd][i] == inf
                                               for mthd in mthds)]
            #gap_sgms = {mthd_res_strs[mthd]: reduce(mul, [gaps[mthd][i]+shift for i in Is], 1)**(1/len(Is))-shift
            #            for mthd in mthds}
            gap_sgms = {mthd_res_strs[mthd]: exp(sum(log(gaps[mthd][i]+shift) for i in Is_gap)/len(Is_gap))-shift
                        for mthd in mthds}

        if boxqp_only and callback_included:
            s_shift = max(1e-4, min(sgaps[mthd][i] for mthd in mthds for i in Is))
            Is_gap = [i for i in Is if not any(isnan(sgaps[mthd][i]) or sgaps[mthd][i] == inf
                                               for mthd in mthds)]

            sgap_sgms = {mthd_res_strs[mthd]: exp(sum(log(sgaps[mthd][i]+s_shift) for i in Is_gap)/len(Is_gap))-s_shift
                         for mthd in mthds}

            #bsols = [float(best_sols[long2short(key)]) for key in mdict[mthds[0]]]
    else:
        BB_cls = {mthd_res_strs[mthd]: 0
                  for mthd in mthds}
        TO_cls = {mthd_res_strs[mthd]: 0
                  for mthd in mthds}
        tave_cls = {mthd_res_strs[mthd]: -1
                    for mthd in mthds}
        nsols_cls = {mthd_res_strs[mthd]: -1
                     for mthd in mthds}
        BB_gaps = {mthd_res_strs[mthd]: -1
                   for mthd in mthds}
        if acopf_only or boxqp_only:
            gap_sgms = {mthd_res_strs[mthd]: -1
                        for mthd in mthds}
        if boxqp_only and callback_included:
            sgap_sgms = {mthd_res_strs[mthd]: -1
                         for mthd in mthds}

        if callback_included:
            nQsols_cls = {mthd_res_strs[mthd]: -1
                          for mthd in mthds}

    BB_gap_strs = {mthd_res_strs[mthd]: f'{100*BB_gaps[mthd_res_strs[mthd]]:.4f}%'
                   for mthd in mthds}

    rdict = {'n': n_cls, 'BB': BB_cls, 'BB_gap': BB_gap_strs, 'TO': TO_cls, 'time': tave_cls, 'nsols': nsols_cls}
    if acopf_only or boxqp_only:
        rdict['gap_sgms'] = gap_sgms

    if boxqp_only and callback_included:
        rdict['sgap_sgms'] = sgap_sgms

    if callback_included:
        rdict['nQsols'] = nQsols_cls

    rstr = rf'{out_path}/res_{runstr}_{cls_str[cls]}.txt'
    with open(rstr, 'wt') as out:
        pprint.pprint(rdict, stream=out)

    M = len(mthds)
    Nc = len(Icls[cls])

    nstr = {mthd: mthd_name_strs[mthd].ljust(8, ' ')
            for mthd in mthds}

    t_str = {mthd: (rf'{tave_cls[mthd_res_strs[mthd]]:.2f}'
                    if cls == SLV else rf'{tave_cls[mthd_res_strs[mthd]]:.1f}'
                    if cls == CNT else '-')
             for mthd in mthds}
    if acopf_only:
        gap_str = {mthd: rf'& {gap_sgms[mthd_res_strs[mthd]]:.3f}\%'
                   for mthd in mthds}
    elif boxqp_only:
        gap_str = {mthd: rf'& {gap_sgms[mthd_res_strs[mthd]]:.2f}\%'
                   for mthd in mthds}
    else:
        gap_str = {mthd: '' for mthd in mthds}

    if boxqp_only and callback_included:
        sgap_str = {mthd: rf'& {sgap_sgms[mthd_res_strs[mthd]]:.3f}\%'
                    for mthd in mthds}
    else:
        sgap_str = {mthd: '' for mthd in mthds}

    BB_str = {mthd: rf'\nicefrac{{{BB_cls[mthd_res_strs[mthd]]}}}{{{Nc}}}'
              for mthd in mthds}

    TO_str = {mthd: (rf'\nicefrac{{{TO_cls[mthd_res_strs[mthd]]}}}{{{Nc}}}'
                     if cls in [CNT,SLV] else '-')
              for mthd in mthds}

    hstr = {mthd: ('        ' + rf'\multirow{{{M}}}{{*}}{{{cls_str[cls]}}}'.ljust(26, ' ')
                   if mthd == mthds[0] else '                                  ')
            for mthd in mthds}

    endstr = (r'\midrule' if cls != USLV else r'\bottomrule')

    tabstr += '\n'.join(rf'{hstr[mthd]}& {nstr[mthd]}& {t_str[mthd]} {gap_str[mthd]} {sgap_str[mthd]}'
                        rf' & {BB_str[mthd]} & {TO_str[mthd]} \\'
                        for mthd in mthds) + f'{endstr}\n'




    rf'''
        \multirow\{{{M}}}\{{*}}\{{{cls_str[cls]}}}   & Bin2  & 6.72 & 2.59\% & \nicefrac\{{5}}\{{14}} & - \\
                                  & Bin3  & 6.45 & 1.72\% & \nicefrac\{{9}}\{{14}} & - \\
                                  & Univariate  & 6.19 & \textbf\{{1.69\%}} & \textbf\{{\nicefrac\{{11}}\{{14}}}} & - \\
                                  & NMDT       & 16.74 & 2.95\% & \nicefrac\{{2}}\{{14}} & - \\
                                  & D-NMDT      & \textbf\{{3.75}} & 2.77\% & \nicefrac\{{5}}\{{14}} & - \\\midrule
        \multirow\{{5}}\{{*}}\{{contested}}& Bin2  & \textbf\{{780.7}} & 12.66\% & \nicefrac\{{0}}\{{2}} & \nicefrac\{{1}}\{{2}} \\
                                  & Bin3  & 899.0 & 12.09\% & \nicefrac\{{1}}\{{2}} & \textbf\{{\nicefrac\{{0}}\{{2}}}} \\
                                  & Univariate   & 1112.7 & 12.09\% & \nicefrac\{{1}}\{{2}} & \nicefrac\{{1}}\{{2}} \\
                                  & NMDT       & 1800.0 & 14.08\% & \nicefrac\{{0}}\{{2}} & \nicefrac\{{2}}\{{2}} \\
                                  & D-NMDT      & 1719.4 & \textbf\{{11.83\%}} & \nicefrac\{{1}}\{{2}} & \nicefrac\{{1}}\{{2}} \\\midrule
        \multirow\{{5}}\{{*}}\{{unsolved}} & Bin2  & - & 31.26\% & \nicefrac\{{1}}\{{43}} & - \\
                                  & Bin3  & - & 30.67\% & \nicefrac\{{11}}\{{43}} & - \\
                                  & Univariate       & - & \textbf\{{30.53\%}} & \textbf\{{\nicefrac\{{29}}\{{43}}}} & - \\
                                  & NMDT       & - & 32.72\% & \nicefrac\{{1}}\{{43}} & - \\
                                  & D-NMDT      & - & 31.86\% & \nicefrac\{{1}}\{{43}} & - \\\bottomrule
    '''

print(tabstr)