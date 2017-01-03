import os
import time
def load_zmip(zmip_result_path):
    import re
    column = []
    file = open(zmip_result_path)
    for line in file:
        line = line.replace('NaN','0')
        line = line.replace('\n','')
        column.append(map(float,re.split(r'\t+', line)))
    file.close()
    return column
def order(zmip_natural):    
    zmip_natural.sort(key=lambda x:x[2], reverse=True)
    
def sincronice_mi(matrix_ref, matrix_evol):
    column = []
    column2 = [] 
    for j in matrix_ref:
        pos  = j[0]     
        pos2 = j[1]   
        #value = j[2]
        column.append(j)
        for e in matrix_evol:
            pos_  = e[0]     
            pos2_ = e[1]
            if (pos_==pos and pos2_==pos2):   
                column2.append(e)
                break    
    return  column, column2
def sincronize_natural_evol_alignments(input_folder,output_folder,pattern_array,reg_init,reg_end_back):
    start_time = time.time()
    print "sincronize_natural_evol_alignments"
    for filename in os.listdir(input_folder):
        if filename.endswith(".cluster") & any(r in filename for r in pattern_array):
            with open(output_folder+"/"+filename,'w') as new_file:
                with open(input_folder+"/"+filename) as old_file:
                    for line in old_file:
                        if('>' in line):
                            new_file.write(line)
                        else:
                            new_file.write(line[reg_init:reg_end_back]+'\n')    
            old_file.close()
            new_file.close()
    print "sincronize_natural_evol_alignments"
    print("--- %s seconds ---" % (time.time() - start_time))

def load_contact_map(contact_map_path):
    import numpy as np    
    with open(contact_map_path) as file:
        l = [map(str,line.split(' ')) for line in file ]
        file.close()
        #test=[[1,2,3,4,5,6,7,8],[1,2,3,4,5,6,7,8],[1,2,3,4,5,6,7,8],[1,2,3,4,5,6,7,8],[1,2,3,4,5,6,7,8],[1,2,3,4,5,6,7,8],[1,2,3,4,5,6,7,8],[1,2,3,4,5,6,7,8]]
        a=np.array(l, order='F')
        return a
        
def sincronize_contact_map(contact_map_path, contact_map_output, reg_init, reg_end_back):
    import numpy as np    
    start_time = time.time()
    print "sincronize_contact_map"
   #with open(contact_map_output,'w') as new_file:
    with open(contact_map_path) as old_file:
        l = [map(str,line.split(' ')) for line in old_file ]
        old_file.close()
        #test=[[1,2,3,4,5,6,7,8],[1,2,3,4,5,6,7,8],[1,2,3,4,5,6,7,8],[1,2,3,4,5,6,7,8],[1,2,3,4,5,6,7,8],[1,2,3,4,5,6,7,8],[1,2,3,4,5,6,7,8],[1,2,3,4,5,6,7,8]]
        a=np.array(l, order='F')
        b=a[reg_init:reg_end_back, reg_init:reg_end_back]
        contacts = np.count_nonzero(b == 'true')
        print "Contacts " + str(contacts)
        np.savetxt(contact_map_output, b, delimiter=' ',fmt="%s")
        numrows = len(b)
        numcols = len(b[0])
        print "ROWS " + str(numrows)
        print "COLUMNS " +str(numcols)
    print "sincronize_contact_map"
    print("--- %s seconds ---" % (time.time() - start_time))
    
    
    