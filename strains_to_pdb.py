import numpy as np 
from collections import defaultdict
import sys

"""
Command-line script made for exchanging mapped strains with corresponding strains provided in csv-file as example output.csv

Script arguments:
input file:     pdb file containing data
strain file:    the file containing the new values which should be mapped to the data-set
outfile:        the filename which the result should be saved to
floatprecision: the number of decimal points to be included from the csv file in the resulting pdb. Anything larger than two still needs to be checked with plotting, as format changes from 00.2f to 0.pf

"""

def change_strain(infile,strainfile,outfile,float_p: int):
    print("received float precision: ",float_p)
    if(float_p==None):
        float_p = 2
    #read header
    f = open(infile,'r')
    header = f.readline()
    f.close()

    #read tail
    termination = False
    with open(infile,'r') as openfileobject:
        for line in openfileobject:
            if "TER" in line:
                tail = line
   
   
    #read number of lines in total
    lines_file = open(infile,"r")
    Counter = 0
    Content = lines_file.read()
    CoList = Content.split("\n")
 
    for i in CoList:
        if i:
            Counter += 1

    num_lines = Counter - 4

    def def_value():
        return "Not Present"
    strain_dict = defaultdict(def_value) #format: key: residue_idx, value: strain


    f_np = np.loadtxt(infile,skiprows=1,dtype=str,max_rows=num_lines) #skip header
    

    B = np.loadtxt(strainfile,delimiter=",",dtype=str,skiprows=1)
    B[:,-1][B[:,-1]==""] = "0.000" #fill empty indices with zero    
    B = B.astype('object')
    B[:,0] = B[:,0].astype('int')
    B[:,1:3] = B[:,1:3].astype('str')
    B[:,-1] = B[:,-1].astype('float32')
    B[:,0]+=1 #align indices
    
    for row in B:
        #print(row[1])
        strain_dict[str(row[0])] = row[-1]

    strain_li = np.zeros((f_np.shape[0],1)).astype("float32")#single column
    for i, entry in enumerate(f_np):
        res_idx = entry[5]
        #print(res_idx)
        strain_li[i] = strain_dict[res_idx]

    f_np[:,-2] = np.array(np.squeeze(strain_li))


    



    def array_to_csv_pdb(outfilename,f,header,tail):
        def write_header_and_tail(filename,header,tail):
            header_num_spaces = 80
            with open(filename, 'a') as file: 
                file.write(tail.ljust(80))
                file.write("ENDMDL".ljust(80)+"\n")
                file.write("END".ljust(80)+"\n")

            with open(filename, 'r+') as fp:
                lines = fp.readlines()     # lines is list of line, each element '...\n'
                lines.insert(0, header.ljust(78))  # you can use any index if you know the line index
                fp.seek(0)                 # file pointer locates at the beginning to write the whole file again
                fp.writelines(lines)       # write whole lists again to the same file

        #maintaining spaces for silly pdb format
        #space_dict = defaultdict(def_value)
        #space_dict["ATOM"] = 0 
        #space_dict["NUM"] = 7 #rjust
        #space_dict["LETTERS"] = 4 #need add two spaces before, then ljust 
        #space_dict["ACID"] = 0
        #space_dict["INDEX"] = 4 #rjust
        #space_dict["X"] = 12 #rjust
        #space_dict["Y"] = 8 #rjust
        #space_dict["Z"] = 8 #rjust
        #space_dict["ONES"] = 6 #rjust
        #space_dict["STRAIN"] = 6 #rjust, format %2.f
        #space_dict["KEY"] = 12 #rjust + 2 spaces after
        

        adjust_list = [7,4,0,4,12,8,8,6,6,12]
        if(int(float_p)>3):
            #adjust_list[-1] -= int(float_p)
            adjust_list[-2] += int(float_p) - 3
            adjust_list[-1] -= int(float_p) - 3




        f_tmp = np.copy(f)
        f_tmp = f_tmp.astype("object")
        
        #print(f_tmp[0,:])
        for i in range(f_tmp.shape[0]):
            #print(f_tmp[i,3]," ",f_tmp[i,4])
            f_tmp[i,3] = f_tmp[i,3] + " " + f_tmp[i,4] #re-add "A" to aminoacid
            if(int(float_p)==2):
                f_tmp[i,-2] = f"{float(f_tmp[i,-2]):05.2f}" #reformat correctly
            else:
                f_tmp[i,-2] = "{:.{prec}f}".format(float(f_tmp[i,-2]),prec=float_p) #reformat correctly
            #f_tmp[i,-2] = f"{float(f_tmp[i,-2]):0{float_p+3}.{float_p}f}" #reformat correctly


        #adjust lengths of substrings
        f_np_tmp = np.hstack((f_tmp[:,:4],f_tmp[:,5:]))
        #f_np_tmp = np.delete(f_np_tmp,-2,axis=1)
        for i in range(f_np_tmp.shape[0]):
            #print(f_np_tmp[i,:])
            for col in range(1,f_np_tmp.shape[1]):
                if(col==2):
                    f_np_tmp[i,col] = "  "+f_np_tmp[i,col].ljust(adjust_list[col-1]," ")
                elif(col==f_np_tmp.shape[1]-1):
                    f_np_tmp[i,col] = f_np_tmp[i,col].rjust(adjust_list[col-1]," ") + "  "
                else:   
                    f_np_tmp[i,col] = f_np_tmp[i,col].rjust(adjust_list[col-1]," ")

       
        np.savetxt(outfilename,f_np_tmp,delimiter="",fmt="%s")
        write_header_and_tail(outfilename,header,tail)
        print("Resulting file written to: ",outfile)


    
        return None
    
    array_to_csv_pdb(outfile,f_np,header,tail)

    return None


if __name__ == "__main__":
    if(sys.argv[1]=="help" or sys.argv[1]=="-h" or sys.argv[1]=="--help" or sys.argv[1]=="-help"):
        print("Calling sequence: \n\t[NECESSARY FLAGS] python strains_to_pdb.py pdb_file_name.pdb csv_file_name.csv outfile_name.pdb \n\t[OPTIONAL FLAGS] decimal precision (integer).\n\tNote: anything else than 2 needs to be checked with plotting procedure, \n\tas float number format changes from 00.2f to 0.(precision)f\n\tDEPENDENCIES: numpy and sys. \n\tYou may install dependencies by issuing: pip install -r requirements.txt")
    else:
        try: 
            float_p = sys.argv[4]
        except IndexError:
            float_p = None
        change_strain(sys.argv[1],sys.argv[2],sys.argv[3],float_p)



