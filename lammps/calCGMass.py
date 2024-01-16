import numpy as np
import sys

#from Bio import SeqIO



def readFasta(fastaFile,start,end):
    #only read single array
    fas=[]
    for seq_record in SeqIO.parse(fastaFile, "fasta"):
        #print(seq_record.id)
        #print(repr(seq_record.seq))
        #print(len(seq_record))
        #fix the start and end issure
        
        for i in range(start-1,end):
            #print(seq_record.seq[i])
            fas.append(seq_record.seq[i])
    
    return fas


def readBoundari(boundaryFile):
    boundry=np.loadtxt(boundaryFile,delimiter=',',dtype=int)
    return boundry


def MassMap():
    aminoacidMass={

        "A":73.03711,
        "R":156.1875,
        "N":114.1038,
        "D":115.0886,
        "C":103.1388,
        "E":129.1155,
        "Q":128.1307,
        "G":57.0519,
        "H":137.1411,
        "I":113.1594,
        "L":113.1594,
        "K":128.1741,
        "M":131.1926,
        "F":147.1766,
        "P":97.1167,
        "S":87.0782,
        "T":101.1051,
        "W":186.2132,
        "Y":163.1760,
        "V":99.1326,
    }

    return aminoacidMass



def CAmass(massMap,fastaSequence, boundary):
    
    MList=[]

    for i in range(len(boundary)):
        mass=0
        if(i==0):
            for j in range(0,boundary[i]+1):
                mass=mass+massMap[fastaSequence[j]]
        elif(i>0):
            for j in range(boundary[i-1]+1,boundary[i]+1):
                mass=mass+massMap[fastaSequence[j]]
        MList.append(mass)
    return MList

def randCAMass(num):
    
    MList=[]

    for i in range(num):
        mass=1
        MList.append(mass)
    return MList



def SaveMass(List,fileName):
    with open(fileName, 'w') as fp:
        for item in List:
            # write each item on a new line
            fp.write("%s\n" % item)
    print('Done')



    
if __name__ == "__main__":
    
    #fastafile name, aminoacid start index(from 1 not 0), aminoacid end index, 
    (junk,num,massFile)=sys.argv
    
    num=int(num)
    #print(MList)
    MList=randCAMass(num=num)

    #save mass into file 
    SaveMass(MList,massFile)
    
    

