############################
## usage: an input file

## aglorithm:
## 0. preparation
## 1. run CG simulations
## 2. analyze CG fluctuation and radial distribution
## 3. compare to AA fluctuation (or radial distribution)
## 4. update CG parameters
## 5. if convergence, finish iteration; else, start a new iteration from step 1.

############################

import sys
import numpy as np
import CG_sim
import random

class params_class:
    def __init__(self, pairstyle):
        self.pairstyle = pairstyle
        self.ps_cutoff = 35
        self.pairs=[]  # list of parameters
        self.npair_types = 0

        self.bondstyle = "harmonic"
        self.bonds=[] # list of parameters
        self.nbond_types = 0

        self.mass=[] # mass list, length = ncg in one molecule
        self.ncg_types = 0  # in fact, length of mass list
        
        
    def init_bond_list(self, infile):
        fh = open(infile)
        for line in fh.readlines():
            if "#" in line:
                continue
            else:
                temp = line.split()
                #bond list format...
                self.bonds.append([ float(temp[3]), float(temp[4])] )
        fh.close()
             
        self.nbond_types = len(self.bonds)

        return True
        

    def init_pair_list(self, infile):
        fh = open(infile)
        for line in fh.readlines():
            if line[0:1]=="#": continue
            temp = line.split()
            if self.pairstyle == "morse":
                self.pairs.append([ float(temp[1]), float(temp[2]), float(temp[3]), float(temp[4]), float(temp[5]), float(temp[6])] )
            else:
                self.pairs.append([ float(temp[1]), float(temp[2]), float(temp[3]), float(temp[4]), float(temp[5]) ] )
        fh.close()
        
        self.npair_types = len(self.pairs)

        return True
        

    def init_mass_list(self, infile):
        fh = open(infile)
        for line in fh.readlines():
            if line[0:1]=="#": continue
            temp = line.split()
            self.mass.append( float(temp[1]) )
        fh.close()
           
        self.ncg_types = len(self.mass)
         
        return True
         
#--------------------------------------------------------
class cg_class:
    
    def __init__(self, index, molid, cgid, cgtype, charge, mass, coord):
        self.index = index
        self.molid = molid
        self.cgid = cgid
        self.cgtype = cgtype
        self.charge = charge
        self.mass = mass
        self.xyz = coord

#--------------------------------------------------------
class topology_class:
    
    def __init__(self, box_list):
        self.natom = 0
        self.nbond =0
        self.nangle = 0
        self.ndihe = 0
        self.nmol = 0
        self.xlim = [box_list[0], box_list[1]]  # xmin, xmax
        self.ylim = [box_list[2], box_list[3]]
        self.zlim = [box_list[4], box_list[5]]
        self.cg_topo = []  # topo for one molecule


    def read_topo(self, infile, massfile):
        
        mass = []
        fmass = open(massfile)
        for line in fmass.readlines():
            mass.append( float(line.split()[1]) )
        fmass.close()
        ncg_per_mol = len(mass)

        index = 0
        molid = 0
        cgid = 0
        ftop = open(infile)
        for line in ftop.readlines()[2:]:  # xyz format
            temp = line.split()
            if cgid < ncg_per_mol:
                coord = np.array( [float(temp[1]), float(temp[2]), float(temp[3]) ] )
                cg = cg_class(index, molid, cgid, cgid, 0.0, mass[cgid], coord)
            else:
                molid = molid + 1
                cgid = 0
                coord = np.array( [float(temp[1]), float(temp[2]), float(temp[3]) ] )
                cg = cg_class(index, molid, cgid, cgid, 0.0, mass[cgid], coord)

            index = index + 1
            cgid = cgid + 1
            self.cg_topo.append(cg)
        ftop.close()
            
        self.natom = index
        self.nmol = molid+1
        self.nbond = int(ncg_per_mol*(ncg_per_mol-1)*0.5*self.nmol)   # cg's within the same mol are all bonded
        
        return True
        
##--------MAIN--------------            
               
def main():
    
    # preparation
    
    debug = False
  
    (junk, cgxyz, blist, nblist, massfile, cgpdb) = sys.argv

    if debug:
        print(cgxyz)
        print(cgpdb)
        print(blist)
        print(nblist)
        print(massfile)
    
    params = params_class("morse")
    params.init_bond_list(blist)
    params.init_pair_list(nblist)
    params.init_mass_list(massfile)
    
    topo = topology_class([-200, 200, -200, 200, -200, 200])
    topo.read_topo(cgxyz, massfile)

    if topo.nbond != len(params.bonds)*topo.nmol:
        print(" Error: mismatch bond list lengths", topo.nbond, len(params.bonds)*topo.nmol)
        print(" Quit ...")
        sys.exit(1)
    elif topo.natom != params.ncg_types*topo.nmol:
        print(" Error: mismatch atom numbers", topo.natom, params.ncg_types*topo.nmol)
        print(" Quit ...")
        sys.exit(1)

    if debug:
        print("Parameters:")
        print(params.pairstyle)
        print(params.ps_cutoff, params.npair_types)
        for ele in params.pairs:
            print(ele)

        print(params.bondstyle)  
        print(params.nbond_types)
        for ele in params.bonds:
            print(ele)

        print(params.ncg_types)
        for ele in params.mass:
            print(ele)
        print
            
        print("Topology:")
        print(topo.natom, topo.nbond, topo.nangle, topo.ndihe, topo.nmol)
        print(topo.xlim, topo.ylim, topo.zlim)
        for ele in topo.cg_topo:
            print(ele.index, ele.molid, ele.cgid, ele.cgtype, ele.charge, ele.mass, ele.xyz)

    infile = "cg_test.in"
    CG_sim.run_cg_sim(infile, topo, params)

    (best_ave, best_std) = CG_ana.get_rmsd()
    print("-1", best_ave, best_std)
    for ele in params.pairs:
        print("%4d%4d%12.4f%12.4f%12.4f"%(ele[0], ele[1], ele[2], ele[3], ele[4]))
    

    # Try to vary parameters for better fit of fluctuations
    max_iter = 0

    for niter in range(max_iter):
             
        # for ele in params.pairs:
        #    print "%4d%4d%12.4f%12.4f%12.4f"%(ele[0], ele[1], ele[2], ele[3], ele[4])

        apair = random.randint(0, len(params.pairs)-1)
        #apair = random.randint(0, 1)
        para2change = random.choice( [ 
            [2, 1.1], [2, 1.08], [2, 1.06], [2, 1.04], [2, 1.03], [2, 1.02], [2, 0.99], [2, 0.97], [2, 0.95], [2, 0.93], [2, 0.9], 
            [3, 1.1], [3, 1.08], [3, 1.06], [3, 1.04], [3, 1.03], [3, 1.02], [3, 0.99], [3, 0.97], [3, 0.95], [3, 0.93], [3, 0.9], 
            [4, 1.1], [4, 1.08], [4, 1.06], [4, 1.04], [4, 1.03], [4, 1.02], [4, 0.99], [4, 0.97], [4, 0.95], [4, 0.93], [4, 0.9] 
            ])

        print("old val:", apair, params.pairs[apair])

        oldval = params.pairs[apair]
        params.pairs[apair][para2change[0]] = params.pairs[apair][para2change[0]] * para2change[1]

        print("new val:", apair, params.pairs[apair])
        
        CG_sim.run_cg_sim(infile, topo, params)
        
        (ave, std) = CG_ana.get_rmsd()

        print(niter, ave, std, apair, para2change)

        if ave < best_ave:
            best_ave = ave
            best_std = std

            fnew = open("new_morselist", "w")
            fnew.write("#    typ1  typ2    D0       alpha       r0              cutoff\n")
            ipt = 1
            for ele in params.pairs:
                print("%4d%4d%12.4f%12.4f%12.4f"%(ele[0], ele[1], ele[2], ele[3], ele[4]))
                fnew.write("%4d%4d%4d%12.4f%12.4f%12.4f%12.4f\n"%(ipt, ele[0], ele[1], ele[2], ele[3], ele[4], ele[5]))
                ipt = ipt + 1
            print("\n")
            fnew.close()
        else:
            params.pairs[apair] = oldval
            print("\n")
    return True

#--------------------------------------

if __name__ == "__main__":
    main()
print("Normal Termination")

