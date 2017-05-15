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
import zipfile
import pandas
import logging
from sklearn import metrics
def getAUC(target,scores):
    fpr, tpr, _ = metrics.roc_curve(target, scores)
    auc = metrics.auc(fpr, tpr)
    return auc
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

'''
Sincroniza y corta los msa evolucionados teniendo en cuenta solo la escructura en comun que existe entre los pdb de la familia
Recibe el pdf recortado con las posiciones que deben mantenerse luego elimina las posiciones del msa evolucionado.
'''
def synchronize_evol_with_cutted_pdb(pdb_complete_path, pdb_cutted_path, clustered_sequences_path, sincronized_evol_path, contact_map_path, sincronized_contact_map):
    start_time = time.time()
    print "sincronize_natural_evol_alignments"  
    df = pandas.read_csv(pdb_complete_path, delim_whitespace=True,header=None)
    df=df.loc[df[4] == 'A']
    df=df.dropna()
    df[5] = df[5].astype('int32')
    df=df.groupby(5).first().reset_index()
    start = df[5].min()
    end = df[5].max()
    df_columnas = pandas.read_csv(pdb_cutted_path, delim_whitespace=True,header=None,usecols=[5])
    df_columnas=df_columnas.dropna()
    df_columnas[5] = df_columnas[5].astype('int32')
    df_columnas=df_columnas.groupby(5).first().reset_index()
    df_columnas[5] = df_columnas[5].apply(lambda x: x - start)
    count=0
    for filename in os.listdir(clustered_sequences_path):
        if filename.endswith(".cluster"):
            with open(sincronized_evol_path+"/"+filename,'w') as new_file:
                with open(clustered_sequences_path+"/"+filename) as old_file:
                    for line in old_file:
                        if('>' in line):
                            line = line.replace('\n','_'+str(count)+'\n')
                            new_file.write(line)
                            count=count+1
                        else:
                            line_array=np.array(list(line))
                            new_line = line_array[df_columnas[5]]
                            new_file.write(new_line.tostring()+'\n')
            old_file.close()
            new_file.close()
    
    
    #Y = np.arange(36).reshape(6,6)
    #test = Y[np.ix_([0,3,5],[0,3,5])]
    cmap = load_contact_map(contact_map_path)
    cmap_sync=cmap[np.ix_(df_columnas[5].tolist(),df_columnas[5].tolist())]
    save_contact_map(cmap_sync, sincronized_contact_map)        
    
    print "sincronize_natural_evol_alignments"
    print("--- %s seconds ---" % (time.time() - start_time))

import random
import string 
def random_char(y):
    return ''.join(random.choice(string.ascii_letters) for x in range(y))

def load_contact_map(contact_map_path, dtype='i4'):
    cmap = np.loadtxt(contact_map_path, dtype=dtype)
    np.set_printoptions(threshold='nan')
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
    
def delete_files(files_pattern):
    files = glob.glob(files_pattern)
    for f in files:
        os.remove(f)
        
def zip_files(files_pattern):
    files = glob.glob(files_pattern)
    for f in files:
        zf = zipfile.ZipFile(f+'.zip', mode='w')
        zf.write(f)       
        
"""
Retorna los PDB/Proteinas a evolucionar.  
Se queda con un pdb por cluster para evitar redundancia.
Con el primer PDB encontrado del cluster
"""        
def find_pdb_to_evolve(family_pdb_information):
    #fields = ["pdb"]
    fields = ['seq',"pdb","chain","cluster","n_residues"]
    df = pandas.read_csv(family_pdb_information,header=0,usecols=fields)
    logging.info("Cantidad Total Proteinas/PDB: " + str(len(df.index)))
    df=df.sort(["cluster","pdb"])
    df=df.groupby("cluster").first()
    df['pdb_folder_name']=df['seq'].str.replace("/","_").str.replace("-","_") + "_" + df['pdb']+"_"+df['chain']
    df['auc']= np.NaN
    df['auc_01']= np.NaN
    df['auc_nat']= np.NaN
    df['auc_nat_01']= np.NaN
    df['spearman_zmip_evol_nat']= np.NaN
    #print df
    logging.info("Cantidad de Proteinas/PDB a evolucionar (Uno por cluster):" + str(len(df.index)))
    print df
    return df
'''
REMOVE PDB HEADER, SAVE ONLY THE ATOMS
'''
def remove_header(pdb_path):
    pdb_temp_path = pdb_path + '_clean'
    with open(pdb_temp_path,'w') as new_file:
        with open(pdb_path) as old_file:
            for line in old_file:
                if(line.startswith("ATOM")):  
                    new_file.write(line)     
    old_file.close()
    new_file.close()
    os.remove(pdb_path)
    os.rename(pdb_temp_path, pdb_path)

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