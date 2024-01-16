import numpy as np
import sys

def readXYZ(frameFile):
    
    # read in all frame coordinates & align each frame to ref
    Nframe = 0  # number of frames
    atom_count = 0
    frame_coord = []  # coordinates for one frame
    frame_list = []   # a list of frames
    type=[]

    frame_fh=open(frameFile)
    
    lines = []
    for line in frame_fh.readlines():
        lines.append(line)
    #print(lines)
    atom_count=int(lines[0])
    print(atom_count)
    Nframe=int(len(lines)/(atom_count+2))
    print(Nframe)
    for i in range(Nframe):
        #print(i)
        for j in range(i*(atom_count+2)+2,(i+1)*(atom_count+2)):
            #print(j)
            temp=lines[j].split()
            frame_coord.append( [ float(temp[1]), float(temp[2]), float(temp[3]) ])
            if(i==0):
                type.append(temp[0])
        frame_list.append(frame_coord)
        frame_coord=[]

    frame_fh.close()
    return frame_list, type




#bulid connection between all CG
def buildAllConnect(xyzFrames,frameIndex,bondCoef,type,nlips):
    #get a frame
    frame=xyzFrames[frameIndex]
    print(type)
    bonds=[]
    
    dis = [[0]*len(frame) for i in range(len(frame))]
    
    for i in range(len(dis)):
        for j in range(len(dis[i])):
            dx=((frame[i][0]-frame[j][0])**2+(frame[i][1]-frame[j][1])**2+(frame[i][2]-frame[j][2])**2)**0.5            
            dis[i][j]=dx
    
    print("done")
    index=1
    
    for i in range(nlips):
        bonds.append([index,bondCoef,dis[(i)*3+1][i*3+2],(i)*3+1,i*3+2])
        index=index+1
        bonds.append([index,bondCoef,dis[(i)*3+2][i*3+3],(i)*3+2,i*3+3])
        index=index+1

    return  bonds,frame

#bulid connection between all CG
def buildAngles(angleCoef,nlips):
    #get a frame
    angles=[]
    print("done")
    index=1

    for i in range(nlips):
        angles.append([index,angleCoef,180,(i)*3+1,i*3+2,i*3+3])
        index=index+1
    
    return  angles


def readMass(mFile):
    Mass=open(mFile)
    
    lines = []
    M=[]
    for line in Mass.readlines():
        lines.append(line)
    
    for i in range(len(lines)):
        temp=lines[i].split()
        #print(temp)
        M.append(float(temp[0]))
    return M


def saveData(fileName, XYZ, bondList,CGMass,anglelist,labellist):


    fo = open(fileName,"w")
    fo.write("# Lammps Data File\n\n")
    fo.write("%8d atoms\n"%len(XYZ))
    #fo.write("%8d atom types\n"%len(XYZ))
    fo.write("%8d atom types\n"%(3))
    fo.write("%8d bonds\n"%len(bondList))
    fo.write("%8d angles\n"%len(anglelist))
    fo.write("\n")

    #bondtype=len(bondList)
    bondtype=1
    fo.write("%8d bond types \n"%bondtype)
    fo.write("\n")

    angletype=1
    fo.write("%8d angle types \n"%angletype)
    fo.write("\n")


    xlo=-200
    ylo=-200
    zlo=-200
    xhi=200
    yhi=200
    zhi=200

    fo.write("%12.3f%12.3f xlo xhi\n"%(xlo, xhi))
    fo.write("%12.3f%12.3f ylo yhi\n"%(ylo, yhi))
    fo.write("%12.3f%12.3f zlo zhi\n"%(zlo, zhi))
    fo.write("\n")
    


    fo.write("Masses\n\n")

    """
    for icg in range(params.ncg_types ):
        fo.write("%8d%12.3f\n"%(icg+1, params.mass[icg]))
    """
    for i in range(3):
    #for i in range(len(CGMass)):
        #fo.write("mass %8d%12.3f\n"%(i+1, CGMass[i]))
        fo.write("%8d%12.3f\n"%(i+1, CGMass[i]))
    fo.write("\n\n")
        

    #fo.write("%8d%12.3f\n"%(1, 1))
    #fo.write("\n")

    #fo.write("bond_style harmonic\n")
    #fo.write("# Bond coeff\n")
    
    fo.write("Bond Coeffs\n\n")
    
    #for ibond in range(params.nbond_types):
    #    fo.write("%8d %12.4f %12.4f\n"%(ibond+1, params.bonds[ibond][0], params.bonds[ibond][1]))
    
    cc=0
    for bb in bondList:
        if(cc==0):
        #fo.write("bond_coeff %8d %12.4f %12.4f\n"%(bb[0],bb[1],bb[2]))
            fo.write("%8d %12.4f %12.4f\n"%(bb[0],bb[1],bb[2]))
        cc=cc+1

    fo.write("\n")

    fo.write("Angle Coeffs\n\n")
    fo.write(("%8d %12.4f %12.4f\n"%(1,10,180)))

    fo.write("\n")

    label={
        "HO1":1,
        "HC1":1,
        "TO1":2,
        "TC1":2,
        "TO2":2,
        "TC2":2,
        "OL":3
    }

    fo.write("Atoms\n\n")
    for i in range(len(XYZ)):
        #format: Index type, molecule-tag atom-type q x y z nx ny nz 
        #(nx,ny,nz are optional - see "true flag" input command
        #index start from 1
        fo.write("%8d"%(i+1))
        moletag=1
        fo.write("%8d"%moletag)
        type=1
        #fo.write("%8d"%(i+1))
        fo.write("%8d"%(label[labellist[i]]))
        charge=0
        fo.write("%12.3f"%charge)

        fo.write("%12.3f"%XYZ[i][0])
        fo.write("%12.3f"%XYZ[i][1])
        fo.write("%12.3f"%XYZ[i][2])
        fo.write("\n")



    fo.write("\n")

    
    fo.write("Bonds\n\n")
    for bb in bondList:
        #fo.write("%8d %8d %8d %8d\n"%(bb[0],bb[0],bb[3],bb[4]))
        fo.write("%8d %8d %8d %8d\n"%(bb[0],1,bb[3],bb[4]))
    fo.write("\n")

##########write angles
    
    fo.write("Angles\n\n")
    for aa in anglelist:
        fo.write("%8d %8d %8d %8d %8d\n"%(aa[0],1,aa[3],aa[4],aa[5])) 
    fo.write("\n")



def SaveInputFile(inpfile,datafileName,bondList):
    fo = open(inpfile,"w")
    fo.write("# test input, created by automated python program\n\n")
    fo.write("units\t real\n")
    fo.write("atom_style\t full\n\n")

    #fo.write("pair_style %s %8.3f\n"%(params.pairstyle, params.ps_cutoff) )

    #non_b='none'
    non_b="lj/charmm/coul/charmm/implicit 8.0 10.0"
    fo.write("pair_style %s\n\n"%non_b)
    
    boundatype='harmonic'
    fo.write("bond_style %s\n\n"%boundatype)

    angletype='harmonic'
    fo.write("angle_style %s\n\n"%angletype)

    fo.write("boundary        p p p\n")

    fo.write("\n")
    #fo.write("region simulation_box block -200 200 -200 200 -200 200")




    #atomtype=len(CGMass)
    bondtype=len(bondList)

    #fo.write("%8d atom types\n"%atomtype)
    #fo.write("bond/types %8d &\n"%bondtype)
    #fo.write("\n")


    fo.write("read_data       "+datafileName+"\n\n")

    fo.write("pair_coeff 1 1 1.00 10.0\n")
    fo.write("pair_coeff 1 2 1.00 10.0\n")
    fo.write("pair_coeff 2 3 1.00 10.0\n")
    fo.write("pair_coeff 2 2 1.00 10.0\n")
    fo.write("pair_coeff 3 3 1.00 10.0\n")
    fo.write("pair_coeff 1 3 1.00 10.0\n\n")

    #for ibond in range(len(params.bonds)):
    #    fo.write("bond_coeff %d %8.4f %8.4f\n"%(ibond+1, params.bonds[ibond][0], params.bonds[ibond][1]))
    #fo.write("\n")

    #fo.write("pair_coeff * * 0.01 0.0001 10.0\n")  # default for pair styles
    #for ipair in range(len(params.pairs)):
    #    fo.write("pair_coeff %d %d %8.4f %8.4f %8.4f\n"%(params.pairs[ipair][0], params.pairs[ipair][1], params.pairs[ipair][2], params.pairs[ipair][3], params.pairs[ipair][4]))
    #fo.write("\n")

    #  fo.write("atom_modify sort 0 0.0\n")  # turn off atom
    """
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
    """    
    fo.write("thermo_style custom step temp ebond epair pe ke etotal press\n")
    fo.write("thermo          100\n")
    fo.write("thermo_modify flush yes\n\n")

    fo.write("dump            1 all dcd 20 cg_test.dcd\n")
    fo.write("dump_modify     1 unwrap yes sort id\n")
    fo.write("dump mydmp all atom 10 dump.min.lammpstrj\n\n")

    fo.write("minimize 10.0 1.0e-4 1000 100000\n")
    fo.write("minimize 0.0 1.0e-8 1000 100000\n")
    fo.write("minimize 1e-12 1.0e-10 1000 10000\n\n")

    fo.write("timestep  1.0\n")
    fo.write("fix       cg_run2  all  nvt  temp 293.0 293.0 100.0\n")
    fo.write("run       100\n")
    fo.write("unfix     cg_run2\n\n")

    fo.write("timestep  2.0\n")
    fo.write("fix       cg_run2  all  nvt  temp 293.0 293.0 100.0\n")
    fo.write("run       150000\n")
    #fo.write("run       1000\n")
    fo.write("unfix     cg_run2\n")
    
    fo.close()




if __name__ == "__main__":
    
    (junk, trjfile,CGmass,datafileName,Inputfile)=sys.argv

    getXYZ,type=readXYZ(trjfile)
    print(getXYZ)
    nlips=3578
    bonds,singleFrame=buildAllConnect(getXYZ,0,10,type,nlips)
    anglelist=buildAngles(10,nlips)


    mass=readMass(CGmass)
    print(mass)
    saveData(datafileName,singleFrame,bonds,mass,anglelist,type)

    SaveInputFile(Inputfile,datafileName,bonds)
