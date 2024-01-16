import os

def savetoXYZ(coordlist,fileName):
    try:
        os.remove(fileName)
    except OSError:
        pass
    
    f=open(fileName,"a")
    f.write(str(len(coordlist))+"\n")
    f.write("txt\n")
    
    for i in range(len(coordlist)):
        f.write(str(coordlist[i]["type"])+"  "+str(coordlist[i]["x"])+"   "+str(coordlist[i]["y"])+"   "+str(coordlist[i]["z"])+"\n")
    
    f.close
