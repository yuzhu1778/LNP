import sys
import os
import time
from subprocess import Popen
#import numpy as np  
#from random import *
#from scipy.optimize import leastsq


#---------------------------------------------------------------------------------

# run lammps CG simulations

def run_cg_sim(inpfile, topo, params, MCstep=40):
        
    #inpfile = "cg_test.in"
    #datfile = "cg_test.dat"
    #logfile = "cg_test.log"

    datfile = inpfile.replace(".in", ".dat")
    logfile = inpfile.replace(".in", ".log")
   
    write_input_file(inpfile, topo, params)
    write_data_file(datfile, topo, params)

    
    proc = Popen(['lmp_mac','-in', inpfile, '-log', logfile, '-screen', 'none'])

    # check if the logfile contains output or not
    normal_termination = False
  
    for i in range(MCstep):
        if normal_termination: 
            break
        else:
            time.sleep(5)
            fh = open(logfile)
            for line in fh.xreadlines():
                if line[0:9] == "Loop time":
                    normal_termination = True
                    break
                elif "ERROR" in line: 
                    os.system('kill -9 %d'%proc.pid)
                    return False
            fh.close()

    #print "the end", normal_termination
    #os.remove(inpfile)
    #os.remove(datfile)
    #os.remove(logfile)
    
    return normal_termination
    
#--------------------------------------------------------------------

def write_input_file(inpfile, topo, params):

    fo = open(inpfile,"w")
    fo.write("# test input, created by automated python program\n\n")
    fo.write("units\t real\n")
    fo.write("atom_style\t full\n\n")

    fo.write("pair_style %s %8.3f\n"%(params.pairstyle, params.ps_cutoff) )
    fo.write("bond_style %s\n\n"%params.bondstyle)

    fo.write("boundary        p p p\n")
    fo.write("read_data       cg_test.dat\n\n")

    #for ibond in range(len(params.bonds)):
    #    fo.write("bond_coeff %d %8.4f %8.4f\n"%(ibond+1, params.bonds[ibond][0], params.bonds[ibond][1]))
    #fo.write("\n")

    fo.write("pair_coeff * * 0.01 0.0001 10.0\n")  # default for pair styles
    for ipair in range(len(params.pairs)):
        fo.write("pair_coeff %d %d %8.4f %8.4f %8.4f\n"%(params.pairs[ipair][0], params.pairs[ipair][1], params.pairs[ipair][2], params.pairs[ipair][3], params.pairs[ipair][4]))
    fo.write("\n")

    #  fo.write("atom_modify sort 0 0.0\n")  # turn off atom
    fo.write("\n")

    for imol in range(topo.nmol):
        molid = imol + 1
        fo.write("group mol%d  molecule %d\n"%(molid, molid))
    fo.write("\n")

    fo.write("neighbor        2.0 bin\n")
    fo.write("neigh_modify    every 1 delay 1 check no\n\n")

    for imol in range(topo.nmol):
        molid = imol + 1
        fo.write("neigh_modify exclude group mol%d mol%d\n"%(molid, molid))
    fo.write("\n")
        
    fo.write("thermo_style custom step temp ebond epair pe ke etotal press\n")
    fo.write("thermo          100\n")
    fo.write("thermo_modify flush yes\n\n")

    fo.write("dump            1 all dcd 20 cg_test.dcd\n")
    fo.write("dump_modify     1 unwrap yes sort id\n\n")

    fo.write("minimize 10.0 1.0e-4 1000 100000\n")
    fo.write("minimize 0.0 1.0e-8 1000 100000\n")
    fo.write("minimize 1e-12 1.0e-10 1000 10000\n\n")

    fo.write("timestep  1.0\n")
    fo.write("fix       cg_run2  all  nvt  temp 293.0 293.0 100.0\n")
    fo.write("run       100\n")
    fo.write("unfix     cg_run2\n\n")

    fo.write("timestep  5.0\n")
    fo.write("fix       cg_run2  all  nvt  temp 293.0 293.0 100.0\n")
    fo.write("run       150000\n")
    #fo.write("run       1000\n")
    fo.write("unfix     cg_run2\n")
    
    fo.close()

    return True

#--------------------------------------------------------------------

def write_data_file(datfile, topo, params):
    
    fo = open(datfile,"w")
    fo.write("# Lammps Data File\n\n")
    fo.write("%8d atoms\n"%topo.natom)
    fo.write("%8d bonds\n"%topo.nbond)
    fo.write("\n")

    fo.write("%8d atom types\n"%params.ncg_types)
    fo.write("%8d bond types\n"%params.nbond_types)
    fo.write("\n")

    fo.write("%12.3f%12.3f xlo xhi\n"%(topo.xlim[0], topo.xlim[1]))
    fo.write("%12.3f%12.3f ylo yhi\n"%(topo.ylim[0], topo.ylim[1]))
    fo.write("%12.3f%12.3f zlo zhi\n"%(topo.zlim[0], topo.zlim[1]))
    fo.write("\n")

    fo.write("Masses\n\n")
    for icg in range(params.ncg_types ):
        fo.write("%8d%12.3f\n"%(icg+1, params.mass[icg]))
    fo.write("\n")

    fo.write("Bond Coeffs\n\n")
    for ibond in range(params.nbond_types):
        fo.write("%8d %12.4f %12.4f\n"%(ibond+1, params.bonds[ibond][0], params.bonds[ibond][1]))
    fo.write("\n")  

    fo.write("Atoms\n\n")
    for icg in range(topo.natom):
        fo.write("%8d"%(icg+1))
        fo.write("%8d"%(topo.cg_topo[icg].molid+1))
        fo.write("%8d"%(topo.cg_topo[icg].cgtype+1))
        fo.write("%12.3f"%topo.cg_topo[icg].charge)
        fo.write("%12.3f"%topo.cg_topo[icg].xyz[0])
        fo.write("%12.3f"%topo.cg_topo[icg].xyz[1])
        fo.write("%12.3f"%topo.cg_topo[icg].xyz[2])
        fo.write("\n")
    fo.write("\n")

    ## need to map a single mol topo to others
    
    fo.write("Bonds\n\n")
    for imol in range(topo.nmol):
        ibond = 0
        for icg in range(params.ncg_types-1):
            for jcg in range(icg+1, params.ncg_types):

                bond_index = imol*params.nbond_types+ibond+1
                bond_type = ibond + 1
                fo.write("%8d %8d %8d %8d\n"%(bond_index, bond_type, imol*params.ncg_types+icg+1, imol*params.ncg_types+jcg+1))
                ibond = ibond + 1
    fo.write("\n")
  
    fo.close()
    return True


