'''
Created on Jan 12, 2017

@author: javi
'''
import numpy as np 
import os
import time
import re
import csv
import glob
def load_zmip(zmip_result_path,window_mi_neightboards=0):
    column = []
    file = open(zmip_result_path)
    for line in file:
        line = line.replace('NaN','0')
        line = line.replace('\n','')
        #column.append(map(float,re.split(r'\t+', line)))
        line_float = map(float,re.split(r'\t+', line))
        if (line_float[0] + window_mi_neightboards >= line_float[1]):
            line_float[2]=0
        column.append(line_float)
    file.close()
    return column
def order(zmip_natural):    
    zmip_natural.sort(key=lambda x:x[2], reverse=True)
    

def save_list_to_csv(filename, zmip, header):
    zmip.insert(0, header)
    with open(filename, "wb") as f:
        writer = csv.writer(f)
        writer.writerows(zmip)    
    zmip.pop(0)
'''
Syncronize matrix_ref with matrix_evol.
'''    
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
def sincronize_natural_evol_msas(input_folder,output_folder,pattern_array,reg_init,reg_end_back):
    start_time = time.time()
    print "sincronize_natural_evol_alignments"
    count=0
    for filename in os.listdir(input_folder):
        if filename.endswith(".cluster") & any(r in filename for r in pattern_array):
            with open(output_folder+"/"+filename,'w') as new_file:
                with open(input_folder+"/"+filename) as old_file:
                    for line in old_file:
                        if('>' in line):
                            line = line.replace('\n','_'+str(count)+'\n')
                            new_file.write(line)
                            count=count+1
                        else:
                            new_file.write(line[reg_init:reg_end_back]+'\n')    
            old_file.close()
            new_file.close()
    print "sincronize_natural_evol_alignments"
    print("--- %s seconds ---" % (time.time() - start_time))

def load_contact_map(contact_map_path):
    cmap = np.loadtxt(contact_map_path, dtype='i4')
    np.set_printoptions(threshold='nan')
    print (cmap)
    return cmap
def load_contact_map_deprecated(contact_map_path):
    with open(contact_map_path) as file:
        l = [map(str,line.split(' ')) for line in file ]
        file.close()
        
        a = np.loadtxt(contact_map_path, dtype='i4')
        np.set_printoptions(threshold='nan')
        #print (cmap)
        #test=[[1,2,3,4,5,6,7,8],[1,2,3,4,5,6,7,8],[1,2,3,4,5,6,7,8],[1,2,3,4,5,6,7,8],[1,2,3,4,5,6,7,8],[1,2,3,4,5,6,7,8],[1,2,3,4,5,6,7,8],[1,2,3,4,5,6,7,8]]
        a=np.array(l, order='F')
        a[a == 'false\n']=0
        a[a == 'false']=0
        a[a == '0\n']=0
        a[a == '0']=0
        a[a == 'true\n']=1
        a[a == 'true']=1
        a[a == '1']=1
        a[a == '1\n']=1
        a[a == '2']=1
        a[a == '2\n']=1
        cmap=np.array(a, dtype='i4')
        #np.set_printoptions(threshold='nan')
        #print (cmap)
        return cmap    
def load_contact_map_(contact_map_path):
    with open(contact_map_path) as file:
        l = [map(str,line.split(' ')) for line in file ]
        file.close()
        #test=[[1,2,3,4,5,6,7,8],[1,2,3,4,5,6,7,8],[1,2,3,4,5,6,7,8],[1,2,3,4,5,6,7,8],[1,2,3,4,5,6,7,8],[1,2,3,4,5,6,7,8],[1,2,3,4,5,6,7,8],[1,2,3,4,5,6,7,8]]
        a=np.array(l, order='F')
        a[a !=0 ] = 1
        cmap=np.array(a, dtype='i4')
        return cmap
        
'''
Sincronize contact map adjusting with the information of reg_init and reg_end_back
contact_map_path: input contact map
contact_map_output: out contact map sincronized
'''        
def sincronize_contact_map(contact_map_path, contact_map_output, reg_init, reg_end_back):
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
    
    
def add_matrix(m1,m2):
    (m, n)=m1.shape
    m3=np.zeros((m, n), dtype='i4')
    for i in range(m):
        for j in range(n):
            m3[i][j] = m1[i][j]+ m2[i][j]
    return m3       

def save_contact_map(m, path):
    np.savetxt(path, m, delimiter=' ',fmt="%s")   
    
def delete_files(folder):
    files = glob.glob(folder+'*')
    for f in files:
        os.remove(f)
'''    
import math, string, sys, fileinput

def range_bytes (): return range(256)
def range_printable(): return (ord(c) for c in string.printable)
def H(data, iterator=range_bytes):
    if not data:
        return 0
    entropy = 0
    for x in iterator():
        p_x = float(data.count(chr(x)))/len(data)
        if p_x > 0:
            entropy += - p_x*math.log(p_x, 2)
    return entropy

def main ():
    for row in fileinput.input():
        string = row.rstrip('\n')
        print ("%s: %f" % (string, H(string, range_printable)))

for str in ['gargleblaster', 'tripleee', 'magnus', 'lkjasdlk',
               'aaaaaaaa', 'sadfasdfasdf', '7&wS/p(']:
    print ("%s: %f" % (str, H(str, range_printable)))
'''    