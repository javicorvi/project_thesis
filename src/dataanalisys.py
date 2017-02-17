'''
Created on Jan 6, 2017

@author: javi
'''
import os
from subprocess import call
import time
import Bio.Cluster
import util
import plot
from sklearn import preprocessing
import numpy as np 
from sklearn.datasets.california_housing import TARGET_FILENAME
 
'''
Calculate the AUC.  
For the protein family (fasta_path) and the contact_map calculates the AUC. 
'''
def auc(fasta_path,contact_map):
    start_time = time.time()
    print "auc"
    #call(["julia", "auc_script.jl" ])
    call(["julia", "mitos/auc.jl",fasta_path,contact_map])
    print "auc"
    print("--- %s seconds ---" % (time.time() - start_time)) 
'''
Calculate the AUC.  
For all the families in clustered_sequences_path and the contact_map_path
The results are:
a zmip for each of the families (result_zmip_path); and
a text file: result_auc_file_name wich contains the AUC for every family. 
'''
def auc_job(pdb_name,model_name,chain_name,contact_map_path,clustered_sequences_path,result_auc_file_name, result_zmip_path):
    start_time = time.time()
    print "auc_process_all"
    #call(["julia", "auc_script.jl" ])
    call(["julia", "mitos/auc_process_all.jl",pdb_name,model_name,chain_name,contact_map_path,clustered_sequences_path,result_auc_file_name,result_zmip_path])
    print "auc_process_all"
    print("--- %s seconds ---" % (time.time() - start_time)) 

'''
Calculates de buslje09 MI corrected.
For the protein family in falsta_path
The result is stored in zmip_result_path
'''     
def buslje09(fasta_path, zmip_result_path):
    start_time = time.time()
    print "buslje09"
    #call(["julia", "auc_script.jl" ])
    call(["julia", "mitos/buslje09.jl",fasta_path, zmip_result_path])
    print "buslje09"
    print("--- %s seconds ---" % (time.time() - start_time))   
        
def buslje09_(input_folder, zmip_result_path,pattern_array):
    start_time = time.time()
    print "buslje09_"
    for filename in os.listdir(input_folder):
        if filename.endswith(".cluster") & any(r in filename for r in pattern_array):
            buslje09(input_folder + filename , zmip_result_path + "zmip" + filename + ".dat")
    print "buslje09_"
    print("--- %s seconds ---" % (time.time() - start_time))

def run_analisys(zmip_natural_result_path, mi_results_path, pattern_array,contact_map_path,outputpath,window):
    #levanto el zmip natural
    zmip_natural = util.load_zmip(zmip_natural_result_path,window)
    util.order(zmip_natural)
    util.save_list_to_csv(zmip_natural_result_path+"_order.csv", zmip_natural, ['Position 1',' Position 2','ZMIP'])
    contact_map=util.load_contact_map(contact_map_path)
    
    for filename in os.listdir(mi_results_path):
        if filename.endswith(".dat") & any(r in filename for r in pattern_array):
            print " Calculation of : " + filename + " with contact map " + contact_map_path
            #levanto el zmip evolucionado 
            zmip_evol = util.load_zmip(mi_results_path + filename,window)
            #sincronizo las senales de coevolucion para calcular spearman rank correlation
            m,m2=util.sincronice_mi(zmip_natural, zmip_evol)
            m_=[row [2] for row in m]
            m2_=[row[2] for row in m2]
            value_spearman = spearman(m_,m2_)
            print " Spearman total " + str(value_spearman)
            
            
            
            
            #MI PLOT
            #contact_map=util.load_contact_map(contact_map_path)
            #plot.contact_map_(contact_map, outputpath)
            #Test load_contact_map
            #v=contact_map[4][2]
            #v=contact_map[4][1]
            #v=contact_map[4][3]
            v=contact_map[80][81]
            m1_norm,m2_norm=normalice_(m_,m2_)
            m_np = np.c_[ np.asarray(m), np.asarray(m1_norm) ]
            m2_np = np.c_[ np.asarray(m2), np.asarray(m2_norm) ]   
            #array for contact mi matrix comparission
            x_nat_t = []
            y_evol_t = []
            x_nat_f = []
            y_evol_f = []
            #scores mi, for roc_curve 
            scores_nat = []
            scores_evol = []
            #1 contact
            #0 no contact
            y_true = []
            for x, y in map(None, m_np, m2_np):
                #se le resta a uno porque los arrays comienzan en la posicion 0
                pos1 = int(x[0]-1)
                pos2 = int(x[1]-1)
                v = contact_map[pos1][pos2]
                #x[0] = posicion 1
                #x[1] = posicion 2
                #x[2] = mi sin normalizar
                #x[3] = mi normalizado
                scores_nat.append(x[3])
                scores_evol.append(y[3])
                if(  v == 0):
                    x_nat_f.append(x[3])
                    y_evol_f.append(y[3])
                    y_true.append(0)
                else:
                    x_nat_t.append(x[3])
                    y_evol_t.append(y[3])
                    y_true.append(1)
            
            #y_true = np.array([0, 1, 1, 1])
            #y_scores = np.array([0.1, 0.4, 0.35, 0.8])
            #y_true = np.array([0, 0, 0, 1])
            #target = np.array([0.1, 0.4, 0.35, 0.8])
           
            labels=['Natural', 'Evol']
            scores=[]
            scores.append(scores_nat)
            scores.append(scores_evol)
            #plot.roc_curve(y_true,scores_nat,scores_evol)
            colors = ['blue', 'red']
            plot.roc_curve(y_true,scores,labels,colors, outputpath+filename+'_roc_curve.png')
            
            plot.contacts_with_mi(x_nat_t,y_evol_t,x_nat_f,y_evol_f,outputpath+filename+'contacts_with_mi.png',filename)
            
            #ordeno zmip evolucionado sincronizado 
            util.order(m2)
            
            util.save_list_to_csv(mi_results_path+filename+"_order.csv", m2, ['Position 1',' Position 2','ZMIP'])
            
            print "TOTAL PAR POSITIONS " + str(len(zmip_natural))
            
            result_file = open(outputpath+filename+".txt","w")
            result_file.write(filename+ '\n')
            result_file.write(" SPEARMAN RANK CORRELATION " + str(value_spearman)+ '\n')
            top_rank(zmip_natural,m2,0.5,contact_map,outputpath+filename+'top_0.5percent_withcon.png',filename,result_file)
            top_rank(zmip_natural,m2,1,contact_map,outputpath+filename+'top_1percent_withcon.png',filename,result_file)
            top_rank(zmip_natural,m2,2,contact_map,outputpath+filename+'top_2percent_withcon.png',filename,result_file)
            top_rank(zmip_natural,m2,3,contact_map,outputpath+filename+'top_3percent_withcon.png',filename,result_file)
            top_rank(zmip_natural,m2,4,contact_map,outputpath+filename+'top_4percent_withcon.png',filename,result_file)
            top_rank(zmip_natural,m2,5,contact_map,outputpath+filename+'top_5percent_withcon.png',filename,result_file)
            result_file.close()
            print '************************************************************************'
            
            
           
'''
Normalize information m,m2
'''
def normalice_(m,m2):
    #import numpy as np   
    #X_train = np.array([[ 1., -1.,  2.],[ 2.,  0.,  0.],[ 0.,  1., -1.]])
    min_max_scaler = preprocessing.MinMaxScaler()
    m_train_minmax = min_max_scaler.fit_transform(m)
    m2_train_minmax = min_max_scaler.fit_transform(m2) 
    return m_train_minmax,m2_train_minmax   
def spearman(x,y):
    return 1 - Bio.Cluster.distancematrix((x,y), dist="s")[1][0]

def matches_coevolved_positions(matrix_ref,matrix_evol):
    data = []
    for j in matrix_ref:
        pos  = j[0]     
        pos2 = j[1]   
        for e in matrix_evol:
            pos_  = e[0]     
            pos2_ = e[1]
            if (pos_==pos and pos2_==pos2):
                #j.append(e[2])
                aux=j[:]
                aux.append(e[2])   
                data.append(aux)
                break    
    return  data
def remove_column(matrix, column):
    return [row[:column] + row[column+1:] for row in matrix]

'''
Plots information about the top_rank.
For example give information about the top_rank matches
Plot a matrix with the contact and the high values (top_rank) of the evolution and the natural msa
'''
def top_rank(x,y,top,contact_map,outputpath,filename,result_file):
    num = len(x)*top/100
    a=x[0:int(num)]
    b=y[0:int(num)]
    print a
    print b
    data=matches_coevolved_positions(a,b)
    a=remove_column(a, 2)
    b=remove_column(b, 2)
    x_nat=[]
    y_nat=[]
    x_evol=[]
    y_evol=[]
    evol_contact_pairs=[]
    nat_contact = 0
    evol_contact = 0
    for x, y in map(None, a, b):
        pos1 = int(x[0]-1)
        pos2 = int(x[1]-1)
        v = contact_map[pos1][pos2]
        if(v == 1):
            nat_contact=nat_contact+1
        pos1 = int(y[0]-1)
        pos2 = int(y[1]-1)
        v = contact_map[pos1][pos2]
        if(v == 1):
            evol_contact=evol_contact+1 
            evol_contact_pairs.append(y)       
        
        #se le resta a uno porque los arrays comienzan en la posicion 0
        x_nat.append(int(x[0]-1))
        y_nat.append(int(x[1]-1))
        x_evol.append(int(y[0]-1))
        y_evol.append(int(y[1]-1))
    
    plot.contact_map_with_top_rank_mi(contact_map,  x_nat, y_nat, x_evol,y_evol,outputpath,filename)
    
    #find information about secondary structure.
    #find information functional information about position.
    
    result_file.write("TOP : "  + str(top) + "% PAR POSITIONS : " + str(num)+ '\n')
    result_file.write("NATURAL CONTACTS QUANTITY : " + str(nat_contact) + " - %"+ str(nat_contact*100/num)+ '\n')
    result_file.write("EVOL CONTACTS QUANTITY : " + str(evol_contact) + " - %"+ str(evol_contact*100/num)+ '\n')
    
    result_file.write("MATCH POSITIONS BETWEEN NAT AND EVOL (NO WINDOW) : " + str(len(data))+ '\n')
    result_file.write("MATCH POSITIONS  : " + str(data) + '\n')
    
    data_contact=[]
    for d in data:
        pos1 = int(d[0]-1)
        pos2 = int(d[1]-1)
        v=contact_map[pos1][pos2]
        if(v==1):
            data_contact.append(d)
    result_file.write("MATCH POSITIONS CONTACTS BETWEEN NAT AND EVOL (NO WINDOW) : " + str(len(data_contact))+ '\n')
    result_file.write("MATCH POSITIONS CONTACTS  : " + str(data_contact) + '\n')
    result_file.write("************************************************************************" + '\n')
    
    
    print "TOP : "  + str(top) + "% PAR POSITIONS : " + str(num)
    print "MATCH POSITIONS BETWEEN NAT AND EVOL (NO WINDOW) : " + str(len(data))
    print "NATURAL CONTACTS QUANTITY : " + str(nat_contact) + " - %"+ str(nat_contact*100/num)
    print "EVOL CONTACTS QUANTITY : " + str(evol_contact) + " - %"+ str(evol_contact*100/num)
    #print data


'''
Plots information about the top_rank.
For example give information about the top_rank matches
Plot a matrix with the contact and the high values (top_rank) of the evolution and the natural msa
'''
def top_rank_desa(x,evol1,evol2,top,contact_map,outputpath,filename,result_file):
    num = len(x)*top/100
    a=x[0:int(num)]
    b=evol1[0:int(num)]
    c=evol2[0:int(num)]
    
    data=matches_coevolved_positions(b,c)
    a=remove_column(a, 2)
    b=remove_column(b, 2)
    c=remove_column(c, 2)
    x_nat=[]
    y_nat=[]
    x_evol1=[]
    y_evol1=[]
    x_evol2=[]
    y_evol2=[]
    evol_contact_pairs=[]
    evol_contact_pairs2=[]
    nat_contact = 0
    evol_contact = 0
    evol_contact2 = 0
    for x, y, f in map(None, a, b, c):
        pos1 = int(x[0]-1)
        pos2 = int(x[1]-1)
        v = contact_map[pos1][pos2]
        if(v == 1):
            nat_contact=nat_contact+1
            
        pos1 = int(y[0]-1)
        pos2 = int(y[1]-1)
        v = contact_map[pos1][pos2]
        if(v == 1):
            evol_contact=evol_contact+1 
            evol_contact_pairs.append(y)
        
        pos1 = int(f[0]-1)
        pos2 = int(f[1]-1)
        v = contact_map[pos1][pos2]
        if(v == 1):
            evol_contact2=evol_contact2+1 
            evol_contact_pairs2.append(y)           
        
        #se le resta a uno porque los arrays comienzan en la posicion 0
        x_nat.append(int(x[0]-1))
        y_nat.append(int(x[1]-1))
        
        x_evol1.append(int(y[0]-1))
        y_evol1.append(int(y[1]-1))
        
        x_evol2.append(int(f[0]-1))
        y_evol2.append(int(f[1]-1))
    
    plot.contact_map_with_top_rank_mi(contact_map,  x_nat, y_nat, x_evol1,y_evol1,x_evol2,y_evol2,outputpath,filename)
    
    #find information about secondary structure.
    #find information functional information about position.
    
    result_file.write("TOP : "  + str(top) + "% PAR POSITIONS : " + str(num)+ '\n')
    result_file.write("NATURAL CONTACTS QUANTITY : " + str(nat_contact) + " - %"+ str(nat_contact*100/num)+ '\n')
    result_file.write("EVOL CONTACTS QUANTITY : " + str(evol_contact) + " - %"+ str(evol_contact*100/num)+ '\n')
    
    result_file.write("MATCH POSITIONS BETWEEN NAT AND EVOL (NO WINDOW) : " + str(len(data))+ '\n')
    result_file.write("MATCH POSITIONS  : " + str(data) + '\n')
    
    data_contact=[]
    for d in data:
        pos1 = int(d[0]-1)
        pos2 = int(d[1]-1)
        v=contact_map[pos1][pos2]
        if(v==1):
            data_contact.append(d)
    result_file.write("MATCH POSITIONS CONTACTS BETWEEN NAT AND EVOL (NO WINDOW) : " + str(len(data_contact))+ '\n')
    result_file.write("MATCH POSITIONS CONTACTS  : " + str(data_contact) + '\n')
    result_file.write("************************************************************************" + '\n')
    
    
    print "TOP : "  + str(top) + "% PAR POSITIONS : " + str(num)
    print "MATCH POSITIONS BETWEEN NAT AND EVOL (NO WINDOW) : " + str(len(data))
    print "NATURAL CONTACTS QUANTITY : " + str(nat_contact) + " - %"+ str(nat_contact*100/num)
    print "EVOL CONTACTS QUANTITY : " + str(evol_contact) + " - %"+ str(evol_contact*100/num)
    #print data    

'''
Not in use.
'''    
def kendall(x,y):
    print 1 - Bio.Cluster.distancematrix((x,y), dist="k")[1][0]

    
        
    
