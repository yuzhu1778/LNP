import numpy as np
import  math
import random


#generate a lipid NP
#nlipids: number of lipids
#lipidLength: number of bead per lipid
#ncavitys: number of cavitys inside lipid
#arealDensity: lipids per 
def LNP(nLipids, lipidLength, ncavitys, center, arealDensity, fatdensity, constants,nfat,pbond):
    if(len(constants)!=2):
        print('bond parameters are not correct')
        
    kbond=constants[0]
    lbond=constants[1]
    #the radius of the LNP
    radius=math.sqrt((nLipids/arealDensity)/(4.0*3.1415926))

    #inside number of lipids
    inner=int(nLipids*radius*radius/((radius+(lipidLength*2-1)*lbond)**2))
    #outer side number of lipids
    outer=int(nLipids-inner)
    coord=[]
    s=3.6/math.sqrt(outer)
    length=0
    dz=2.0/outer
    z=1.0-dz/2.0

    #get the position of out layer
    #go through lipids
    for i in range(outer):
        r=math.sqrt(1.0-z*z)

        #go through each lipid
        for j in range(lipidLength):
            #print(j)
            
            p={
            "Record_Type":  "ATOM",
            "Res_Sequence_N":   "0",
            "type": "T",
            "x":    0,
            "y":    0,
            "z":    0,
            "Res_Name": "GLP"
            }
            p["Res_Sequence_N"]=i+1
            p["x"]=center["x"]+r*math.cos(length)*(radius+lbond*j)
            p["y"]=center["y"]+r*math.sin(length)*(radius+lbond*j)
            p["z"]=center["z"]+z*(radius+lbond*j)
            if(j==lipidLength-1):
                #out layer head
                p["type"]="HO1"
            else:
                #out layer tail
                 
                p["type"]="TO"+str(lipidLength-j-1)
            
            coord.append(p)
        z=z-dz
        length=length+s/r

    #create cavities
    #determine number of lipid per cavity
    lip_cavity=[]
    center_cavity=[]
    radius_cavity=[]
    for i in range(ncavitys):
        lip_cavity.append(int(inner/ncavitys))
        
    lip_cavity[-1]=inner-int(inner/ncavitys)*(ncavitys-1)
    
    #calculate the radius of each cavity
    for i in range(ncavitys):
        radius_cavity.append(math.sqrt((lip_cavity[i]/arealDensity)/(4.0*3.1415926)))

    
    #determine the center of each cavity
    
    count=0
    ncavitys=int(0.5*ncavitys)
    while count<ncavitys:
        x=random.uniform(0,1)*(radius-radius_cavity[i]-(lipidLength-1)*lbond)*2-radius+radius_cavity[i]+(lipidLength-1)*lbond
        y=random.uniform(0,1)*(radius-radius_cavity[i]-(lipidLength-1)*lbond)*2-radius+radius_cavity[i]+(lipidLength-1)*lbond
        z=random.uniform(0,1)*(radius-radius_cavity[i]-(lipidLength-1)*lbond)*2-radius+radius_cavity[i]+(lipidLength-1)*lbond

        accept=True
        if math.sqrt(x*x+y*y+z*z)>(radius-radius_cavity[i]-(lipidLength-1)*lbond):
            accept=False

        if(len(center_cavity)==0):
            pass
            #center_cavity.append([x,y,z])
        else:
            for pp in center_cavity:
                dx=pp[0]-x;
                dy=pp[1]-y;
                dz=pp[2]-z;
                
                d=math.sqrt(dx*dx+dy*dy+dz*dz)
                if d<2*radius_cavity[0]+(lipidLength-1)*lbond*2:
                    accept=False
        
        if accept:
            center_cavity.append([x,y,z])
            count=count+1
    
    #print(center_cavity)
    #print(radius_cavity)


    #print(center_cavity)
    
    for i in range(ncavitys):
        s=3.6/math.sqrt(lip_cavity[i])
        length=0
        dz=2.0/lip_cavity[i]
        z=1.0-dz/2.0

        for j in range(lip_cavity[i]):
            r=math.sqrt(1.0-z*z)
            for k in range(lipidLength):
                p={
                "Record_Type":  "ATOM",
                "Res_Sequence_N":   "0",
                "type": "T",
                "x":    0,
                "y":    0,
                "z":    0,
                "Res_Name": "GCA"
                }
                p["Res_Sequence_N"]=outer+i*j+j+1

                p["x"]=center_cavity[i][0]+r*math.cos(length)*(radius_cavity[i]+lbond*k)
                p["y"]=center_cavity[i][1]+r*math.sin(length)*(radius_cavity[i]+lbond*k)
                p["z"]=center_cavity[i][2]+z*(radius_cavity[i]+lbond*k)
                if(k==0):
                    #cavity head
                    p["type"]="HC1"
                else:
                    #cavity tail
                    p["type"]="TC"+str(k)
                coord.append(p)
            z=z-dz
            length=length+s/r

    totalInner=0
    for i in lip_cavity:
        totalInner=totalInner+i
    
    print(outer+totalInner)

    #fill the space with oil particle
    count=0
    while count<nfat:
        p={                
                "Record_Type":  "ATOM",
                "Res_Sequence_N":   "0",
                "type": "T",
                "x":    0,
                "y":    0,
                "z":    0,
                "Res_Name": "OL"
                }
        p["Res_Sequence_N"]=outer+totalInner+count+1
    
        p["x"]=random.uniform(0,1)*(radius)*2-radius
        p["y"]=random.uniform(0,1)*(radius)*2-radius
        p["z"]=random.uniform(0,1)*(radius)*2-radius
        p["type"]="OL"
        
        accept=True
        if math.sqrt(p["x"]*p["x"]+p["y"]*p["y"]+p["z"]*p["z"])>(radius):
            accept=False

        for cen in center_cavity:
                dx=p["x"]-cen[0];
                dy=p["y"]-cen[1];
                dz=p["z"]-cen[2];
                d=math.sqrt(dx*dx+dy*dy+dz*dz)
                if d < radius_cavity[0]+lipidLength*pbond:
                    accept=False
        if accept:  
            coord.append(p)
            count=count+1


    return radius,coord

    