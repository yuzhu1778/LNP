import LNP
import saveCoord

import savepdb

if __name__=='__main__':

    COM={
        "x":0,
        "y":0,
        "z":0
    }

    LNPcons=[100,0.7]

    #since the space is not enough to have enough cavitys, the number of cavitys is half
    r,coor=LNP.LNP(nLipids=5000,lipidLength=3,ncavitys=2,center=COM,arealDensity=3.45,fatdensity=3.45,constants=LNPcons,nfat=4000,pbond=0.7)
    print(r)

    saveCoord.savetoXYZ(coor,"test.xyz")
    
    size=[20*r,20*r,20*r]
    savepdb.savetoPDB(coordlist=coor, fileName="stest.pdb", size=size)