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
TNMDT = Model_qcqp.TNMDT
TDNMDT = Model_qcqp.TDNMDT
BIN2 = Model_qcqp.BIN2
BIN3 = Model_qcqp.BIN3
MCCORMICK = Model_qcqp.MCCORMICK
UNIVAR_NOSTUB = Model_qcqp.UNIVAR_NOSTUB
UNIVAR_NOST = Model_qcqp.UNIVAR_NOST
DIRECT = 9001

### Primary settings; user input. Includes file paths and model depths. ###
# File Paths
in_path = r'./All_instances/original'
out_path = r'./All_instances/testcases_uni_ref'

# Model Depths
L = 3
L1 = 10

# Flags
remodel_anyway = 1  # Re-create any models that already exist.

# List of methods #
#mthds = [BIN2, BIN3, UNIVAR, TNMDT, TDNMDT, NMDT, DNMDT, UNIVAR_NOSTUB]
mthds = [DIRECT, UNIVAR_NOST, UNIVAR_NOSTUB]

### Other model settings; user input ###
# DNMDT #
# ldas: list of lambda values to use for DNMT.
ldas = [0, 1]

# NMDT #
disc_left = 1
disc_right = 0

# Sawtooth #
st_do_approx = 0

### Other setup ###
### No input from users past this point ###
mdict = {mthd: {} for mthd in mthds}
ldas_str = f'ldas={"-".join(str(i) for i in ldas)}'

# NMTT: NMDT, with pure quadratic terms tightened with sawtooth LB improvement
dnmdt_str = f'_DNMDT_L={L}_ldas={ldas_str}'
nmdt_str = f'_NMDT_L={L}_discleft={disc_left}_discright={disc_right}'
tdnmdt_str = f'_T-DNMDT_L={L}_ldas={ldas_str}'
tnmdt_str = f'_T-NMDT_L={L}_discleft={disc_left}_discright={disc_right}'
univar_str = f'_UNIVAR_L={L}_L1={L1}_doapprox={st_do_approx}'
univar_nostub_str = f'_UNIVAR_NOSTUB_L1={L1}'
univar_nost_str = f'_UNIVAR_NOST={L1}'
bin2_str = f'_BIN2_L={L}_L1={L1}'
bin3_str = f'_BIN3_L={L}_L1={L1}'
mc_str = f'_MCCORMICK'
direct_str = f''


mthd_res_strs = {UNIVAR: univar_str[1:],
                 DNMDT: dnmdt_str[1:],
                 NMDT: nmdt_str[1:],
                 TDNMDT: tdnmdt_str[1:],
                 TNMDT: tnmdt_str[1:],
                 BIN2: bin2_str[1:],
                 BIN3: bin3_str[1:],
                 UNIVAR_NOST: univar_nost_str[1:],
                 UNIVAR_NOSTUB: univar_nost_str[1:],
                 MCCORMICK: 'MCCORMICK',
                 DIRECT: 'DIRECT'}

mthd_strs = {UNIVAR: univar_str, DNMDT: dnmdt_str, NMDT: nmdt_str,
             DIRECT: direct_str, TDNMDT: tdnmdt_str, TNMDT: tnmdt_str,
             BIN2: bin2_str, BIN3: bin3_str, MCCORMICK: mc_str,
             UNIVAR_NOSTUB: univar_nostub_str, UNIVAR_NOST: univar_nost_str}

if not (disc_right or disc_left):
    disc_left = 0

# File path
myfpath = rf'{in_path}/*.lp'
fpaths = glb(myfpath)

# Write all the models. Unless remodel_anyway=True, ignore models that already exist
for fpath in fpaths:
    (fp, fnameext) = path.split(fpath)
    (pname, fext) = fsplit(fnameext)
    # fp: original file path
    # fname: name of file, without extension
    # fext: read-in extension

    mpath_exists = 1
    # Write all models for the current problem to files
    dir_str = mthd_strs[DIRECT]
    for mthd in mthds:
        mthd_str = mthd_strs[mthd]
        mpath = rf'{out_path}/{pname}{mthd_str}.lp'
        rpath = rf'{out_path}/Results/{pname}{mthd_str}'

        m = read(fpath)
        m.update()

        if not (path.exists(mpath)) or remodel_anyway:
            if mthd != DIRECT:
                    myL = (2*L if mthd in [NMDT, TNMDT] else L)
                    mdl_holder = Model_qcqp(m, L=myL, L1=L1, ldas=[0], method=mthd,
                                            disc_right=disc_right, disc_left=disc_left, st_do_approx=st_do_approx)

                    mdl = mdl_holder.mdl
                    mdl.update()
            else:
                mdl = m

            # write model to file
            mdl.write(rf'{out_path}/{pname}{mthd_str}.lp')