from glob import glob as glb
from os import path
from gurobipy import *
from Model_qcqp import Model_qcqp

fsplit = path.splitext

# Flags from Model_qcqp
UNIVAR = Model_qcqp.UNIVAR
NMDT = Model_qcqp.NMDT
DNMDT = Model_qcqp.DNMDT
DIRECT = 9001

### File settings; user input ###
p_classes = ['basic', 'small', 'extended', 'extended2']
#p_classes = ['basic']
out_path = r'./boxQP/mdl-data'

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
#mthds = [UNIVAR]

### Other setup ###
ldas_str = f'ldas={"-".join(str(i) for i in ldas)}'
if not (disc_right or disc_left):
    disc_left = 0
for pcls in p_classes:
    # file path per pcls; user input
    fpath = rf'./boxQP/LP-data/{pcls}-LP/*.lp'

    ### No input from users past this point ###
    fpaths = glb(fpath)

    for fpath in fpaths:
        (fp, fnameext) = path.split(fpath)
        (fname, fext) = fsplit(fnameext)
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

        # Write all models for the current problem to files
        for mthd in mthds:
            if mthd != DIRECT:
                myL = (2*L if mthd == NMDT else L)
                mdl_holder = Model_qcqp(m, L=myL, L1=L1, ldas=[0], method=mthd,
                                        disc_right=disc_right, disc_left=disc_left, st_do_approx=st_do_approx)

                mdl = mdl_holder.mdl
                mdl.update()
                if mthd == DNMDT:
                    mdl.write(rf'{out_path}/{pcls}/{fname}_DNMDT_L={L}_ldas={ldas_str}.lp')
                elif mthd == NMDT:
                    mdl.write(rf'{out_path}/{pcls}/{fname}_NMDT_L={L}_discleft={disc_left}_discright={disc_right}.lp')
                elif mthd == UNIVAR:
                    mdl.write(rf'{out_path}/{pcls}/{fname}_UNIVAR_L={L}_L1={L1}_doapprox={st_do_approx}.lp')
                else:
                    continue  # invalid method, so skip it
            else:
                m.write(rf'{out_path}/{pcls}/{fname}.lp')