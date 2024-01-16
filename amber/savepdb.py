import os
import hybrid_36

#unit is angstrom

def savetoPDB(coordlist,size,fileName):
    try:
        os.remove(fileName)
    except OSError:
        pass
    
    f=open(fileName,"a")
    

    f.write("CRYST1    1.000    1.000    1.000  "+str(size[0])+"  "+str(size[1])+"  "+str(size[2])+" P 1\n")
    
    #for i in range(1):
    for i in range(len(coordlist)):
        j={
            "Record_Type": coordlist[i]["Record_Type"],
            "Atom_Serial_Number":   hybrid_36.hy36encode(5,i+1),
            "Atom_Name":    coordlist[i]["type"],
            "Res_Name":     coordlist[i]["Res_Name"],
            "Res_Number":   hybrid_36.hy36encode(4,coordlist[i]["Res_Sequence_N"]),
            "X":    str(coordlist[i]["x"]*10.0),
            "Y":    str(coordlist[i]["y"]*10.0),
            "Z":    str(coordlist[i]["z"]*10.0),
            "Occupancy":    "1.00",
            "Temperature_factor":   "0.00",
            "Chain_Identifier":     "A",
            "Element_Symbol":   "CG",
        }

        j["Record_Type"]=j["Record_Type"].ljust(6)
        j["Atom_Serial_Number"]=j["Atom_Serial_Number"].rjust(5)
        j["Atom_Name"]=j["Atom_Name"].ljust(4)
        j["Res_Name"]=j["Res_Name"].ljust(3)
        j["Res_Number"]=j["Res_Number"].rjust(4)

        j["X"]=str('%8.3f' % (float(j["X"]))).rjust(8)
        j["Y"]=str('%8.3f' % (float(j["Y"]))).rjust(8)
        j["Z"]=str('%8.3f' % (float(j["Z"]))).rjust(8)
        
        j["Occupancy"]=j["Occupancy"].rjust(6)
        j["Temperature_factor"]=j["Temperature_factor"].rjust(6)
        j["Element_Symbol"]=j["Element_Symbol"].rjust(2)
        #f.write(j[0]+j[1]+" "+j[2]+"     "+j[4]+j[3]+"    "+j[6]+j[7]+j[8]+"\n")
        
        s=(j["Record_Type"]+
            j["Atom_Serial_Number"]+" "+
            j["Atom_Name"]+" "+
            j["Res_Name"]+" "+
            j["Chain_Identifier"]+
            j["Res_Number"]+
            "    "+j["X"]+j["Y"]+j["Z"]+
            j["Occupancy"]+
            j["Temperature_factor"]+
            "          "+j["Element_Symbol"]+
            "\n")


        f.write(s)
    f.write("END")

    f.close
