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
import glob
import msa
from sklearn.datasets.california_housing import TARGET_FILENAME
import constants as cons
import pandas 
import logging
'''
Calculate the AUC.  
For the protein family (fasta_path) and the contact_map calculates the AUC. 
'''
#interpolate false
def auc(fasta_path,contact_map):
    start_time = time.time()
    print "auc"
    #call(["julia", "auc_script.jl" ])
    call([cons.julia_exe_path, "mitos/auc.jl",fasta_path,contact_map])
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
    call([cons.julia_exe_path, "mitos/auc_process_all.jl",pdb_name,model_name,chain_name,contact_map_path,clustered_sequences_path,result_auc_file_name,result_zmip_path])
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
    #call(["julia"])
    #call(["julia04/bin/julia", "mitos/buslje09.jl" ])
    #h=heatmap(full(zmip), yflip=true)
    #call(["julia04/bin/julia", "mitos/buslje09.jl",fasta_path, zmip_result_path])
    call([cons.julia_exe_path, "mitos/buslje09.jl",fasta_path, zmip_result_path])
    print "buslje09"
    print("--- %s seconds ---" % (time.time() - start_time))   
        
def buslje09_(input_folder, zmip_result_path,pattern_array=["sequences"]): 
    start_time = time.time()
    print "buslje09_"
    for filename in os.listdir(input_folder):
        if filename.endswith(".cluster") & any(r in filename for r in pattern_array):
            buslje09(input_folder + filename , zmip_result_path + "zmip" + filename + ".dat")
    print "buslje09_"
    print("--- %s seconds ---" % (time.time() - start_time))

def run_analisys(df,index, zmip_natural_result_path, mi_results_path, pattern_array,contact_map_path,outputpath,window):
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
            df.set_value(index, 'spearman_zmip_evol_nat', value_spearman) 
            #MI PLOT
            #contact_map=util.load_contact_map(contact_map_path)
            #plot.contact_map_(contact_map, outputpath)
            #Test load_contact_map
            #v=contact_map[4][2]
            #v=contact_map[4][1]
            #v=contact_map[4][3]
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
            plot.roc_curve(df,index,y_true,scores,labels,colors, outputpath+filename+'_roc_curve.png')
            
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

def run_analisys_singular(df,index, zmip_natural_result_path, mi_result_file_path, pattern_array,contact_map_path,outputpath,window):
    #levanto el zmip natural
    zmip_natural = util.load_zmip(zmip_natural_result_path,window)
    util.order(zmip_natural)
    util.save_list_to_csv(zmip_natural_result_path+"_order.csv", zmip_natural, ['Position 1',' Position 2','ZMIP'])
    contact_map=util.load_contact_map(contact_map_path)
    
    print " Calculation of : " + filename + " with contact map " + contact_map_path
    #levanto el zmip evolucionado 
    zmip_evol = util.load_zmip(mi_results_path + filename,window)
    #sincronizo las senales de coevolucion para calcular spearman rank correlation
    m,m2=util.sincronice_mi(zmip_natural, zmip_evol)
    m_=[row [2] for row in m]
    m2_=[row[2] for row in m2]
    value_spearman = spearman(m_,m2_)
    df.set_value(index, 'spearman_zmip_evol_nat', value_spearman) 
    #MI PLOT
    #contact_map=util.load_contact_map(contact_map_path)
    #plot.contact_map_(contact_map, outputpath)
    #Test load_contact_map
    #v=contact_map[4][2]
    #v=contact_map[4][1]
    #v=contact_map[4][3]
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
            
        labels=['Natural', 'Evol']
        scores=[]
        scores.append(scores_nat)
        scores.append(scores_evol)
        #plot.roc_curve(y_true,scores_nat,scores_evol)
        colors = ['blue', 'red']
        plot.roc_curve(df,index,y_true,scores,labels,colors, outputpath+filename+'_roc_curve.png')
            
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
def normalice_desarrollo(m,m2,m3):
    #import numpy as np   
    #X_train = np.array([[ 1., -1.,  2.],[ 2.,  0.,  0.],[ 0.,  1., -1.]])
    min_max_scaler = preprocessing.MinMaxScaler()
    m_train_minmax = min_max_scaler.fit_transform(m)
    m2_train_minmax = min_max_scaler.fit_transform(m2) 
    m3_train_minmax = min_max_scaler.fit_transform(m3) 
    return m_train_minmax,m2_train_minmax,m3_train_minmax           
           
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
    
    plot.contact_map_with_top_rank_mi_desarrollo(contact_map,  x_nat, y_nat, x_evol1,y_evol1,x_evol2,y_evol2,outputpath,filename)
    
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
Toma todas las matrices de contactos de todas la proteinas evolucionadas,  retorna una matriz con la sumas y otra matriz con las probabilidad de los contactos.
'''
def sum_contact_map(family_folder,pdb_to_evol_df):
    family_folder_pdb = family_folder+"/PDB/"
    cmap_sum = None
    for index,pdb_protein_to_evolve in pdb_to_evol_df.iterrows():
        pdb_folder = family_folder_pdb + pdb_protein_to_evolve['pdb_folder_name']
        if(os.path.isdir(pdb_folder)): 
            contact_map = pdb_folder + "/contact_map.dat"
            cmap = util.load_contact_map(contact_map)
            print contact_map
            if(cmap_sum==None):
                cmap_sum = cmap
            else:
                if(cmap.shape==cmap_sum.shape):
                    cmap_sum = cmap_sum + cmap
                else:
                    print " diferent size natural : " + str(cmap_sum.shape)  + " AND " + contact_map + " : " + str(cmap.shape)
                    #pdb_to_evol_df=pdb_to_evol_df.drop(index)
    #print cmap_sum 
    util.save_contact_map(cmap_sum, family_folder + "/sum_contact_map.dat")
    cmap_sum = cmap_sum.astype(float)
    camp_prob = cmap_sum / pdb_to_evol_df.shape[0]
    #print camp_prob   
    util.save_contact_map(camp_prob, family_folder + "/prob_contact_map.dat") 
    plot.contact_map(camp_prob,family_folder + "/prob_contact_map.png")
    conserved_contacts = np.count_nonzero(camp_prob == 1.0)  
    print conserved_contacts
"""
Lee la informacion sobre consevacion (KL) por columna de cada una de las proteinas evolucionadas. 
No se esta aplicando ningun entrecruzamiento de la informacion.
Solamente se esta ploteando la conservacion por columna para cada una de las proteinas (el grafico no puede apreciar resultados concretos)  
"""    
def comparative_conservation(family_folder, family_name, pdb_to_evol_df):
    natural_msa_conservation= family_folder + "/"+family_name+".fasta_data_kl.csv"
    family_folder_pdb = family_folder+"/PDB/"
    #[ name for name in os.listdir(thedir) if os.path.isdir(os.path.join(thedir, name)) ]
    msas_entropy=[]
    msa_entropy_media=[]
    df=msa.read_conservation(natural_msa_conservation)
    df = df.dropna()
    msa_entropy = [df['Entropy'].tolist(),family_name + "_NATURAL"]
    msas_entropy.append(msa_entropy)
    cant=0
    for index,pdb_protein_to_evolve in pdb_to_evol_df.iterrows():
        pdb_folder = family_folder_pdb + pdb_protein_to_evolve['pdb_folder_name']
        if(os.path.isdir(pdb_folder)): 
            conservation_file = pdb_folder + "/clustered_sequences/information/*kl.csv"
            conservation_file = glob.glob(conservation_file) 
            df=msa.read_conservation(conservation_file[0])
            df = df.dropna()
            msa_entropy = [df['Entropy'].tolist(),pdb_folder]
            msas_entropy.append(msa_entropy)
            if(cant==0):
                msa_entropy_media = df['Entropy'].tolist()
            else:
                msa_entropy_media = [x + y for x, y in zip(msa_entropy_media , df['Entropy'].tolist() )]
            cant=cant+1
    msa_entropy_media =   [x / cant  for x in msa_entropy_media]     
    plot.conservation_between_msas(msas_entropy,family_folder + "/conservation.png")  
    #Media Graphic with natural
    msas_entropy=[]
    df=msa.read_conservation(natural_msa_conservation)
    df = df.dropna()
    msa_entropy = [df['Entropy'].tolist(),family_name + "_NATURAL"]
    msas_entropy.append(msa_entropy)
    msas_entropy.append([msa_entropy_media,"MEDIA"])
    plot.conservation_between_msas(msas_entropy,family_folder + "/conservation_media.png") 
"""
Esta funcion toma el top de MI de todas las proteinas evolucionadas y luego realiza una agrupacion indicando la cantidad de veces que aparecen los pares.
Ordena los pares de forma descendente, osea los pares que mas aparecen en el top quedan arriba. 
Ademas se agrega la columna indicando la probabilidad de contacto que existen entre ellos.
"""
def comparative_mi_information(family_folder, family_name,top, window, pdb_to_evol_df):     
    import matplotlib.pyplot as plt      
    logging.info('Begin of the execution process family MI information')
    family_folder_pdb = family_folder+"/PDB/"
    #for protein_pdb in os.listdir(family_folder_pdb):
    fields=["Position1","Position2","Count"]
    df_total = pandas.DataFrame([],columns=fields)
    cant = 0
    for index,pdb_protein_to_evolve in pdb_to_evol_df.iterrows():
        pdb_folder = family_folder_pdb + pdb_protein_to_evolve['pdb_folder_name']
        if(os.path.isdir(pdb_folder)): 
            zmip_file_pattern = pdb_folder + "/mi_data/zmipsequences*.dat"
            zmip_file = glob.glob(zmip_file_pattern) 
            if(len(zmip_file)==0):
                logging.error('No existe archivo de informacion mutua de la familia ' + pdb_folder)
                return
            if(len(zmip_file)>1):
                logging.error('Existe mas de un archivo de informacion mutua de la familia ' + pdb_folder)
                return
            cant = cant + 1
            zmip_evol = util.load_zmip(zmip_file[0],window)
            util.order(zmip_evol)
            num = len(zmip_evol)*top/100
            zmip_evol_top=zmip_evol[0:int(num)]
            df = pandas.DataFrame(zmip_evol_top,columns=fields)
            df_total=df_total.append(df)
            #print df_total
    print df_total
    
    """a=[14.0,76.0, 45.2345]
    b=[23.0,54.0, 34.5]
    c=[[14.0,76.0, 45.2345],[23.0,54.0, 34.5],[45.0,90.0, 34.5]]
    #c.append(a)
    #c.append(b)
    df = pandas.DataFrame(c,columns=fields)
    d=[[11.0,76.0, 45.2345],[23.0,54.0, 34.5],[45.0,90.0, 34.5]]
    #d.append(a)
    #d.append(b)
    df2 = pandas.DataFrame(d,columns=fields)
    dfres=df.append(df2)
    e=[[10.0,76.0, 45.2345],[23.0,54.0, 34.5],[45.0,90.0, 34.5]]
    dfe = pandas.DataFrame(e,columns=fields)
    dfres=dfres.append(dfe)
   
    """
    #After append all the evolution MI TOP do this
    #counts = dfres.groupby(['Position1','Position2']).size()
    #print counts
    counts_df = pandas.DataFrame(df_total.groupby(['Position1','Position2']).size().rename('Count'))
    print counts_df
    sorted_df=counts_df.sort_values(by=['Count'],ascending=[False])
    print sorted_df
    sorted_df['ProbContact']=pandas.Series(0.0, index=sorted_df.index)
    sorted_df['ProbTop']=pandas.Series(0.0, index=sorted_df.index)
    prob_contact_map = util.load_contact_map(family_folder + "/prob_contact_map.dat",np.float64)
    print prob_contact_map
    for index,mi_par in sorted_df.iterrows():
        #por bug arreglar
        #print mi_par
        pos1 = int(index[0]-1)
        pos2 = int(index[1]-1)
        v=prob_contact_map[pos1][pos2]
        sorted_df.set_value(index, 'ProbContact' , v)
        prob_top = mi_par['Count'] * 100 / cant 
        sorted_df.set_value(index, 'ProbTop' , prob_top/100)
        #sorted_df[index]['ProbContact']=v
    
    
    sorted_df.to_csv(family_folder + "/top_family_mi.csv", sep='\t', encoding='utf-8')
    
    
    correlation_p = sorted_df['ProbContact'].corr(sorted_df['ProbTop'], method='pearson')
    correlation_k = sorted_df['ProbContact'].corr(sorted_df['ProbTop'], method='kendall')
    correlation_s = sorted_df['ProbContact'].corr(sorted_df['ProbTop'], method='spearman')
    
    '''correlation_sp = df.corr(method='pearson')
    correlation_sp = df.corr(method='kendall')
    correlation_sp = df.corr(method='spearman')
    '''
    mean = sorted_df["ProbTop"].mean()
    median = sorted_df["ProbTop"].median()
    var = sorted_df["ProbTop"].var()
    mode = sorted_df["ProbTop"].mode()
    sorted_df.plot.scatter(x='ProbTop', y='ProbContact');
    plt.savefig(family_folder + "/top_family_mi.png");
    plt.show()
    plt.gcf().clear()
    print sorted_df
    

                     
'''
Not in use.
'''    
def kendall(x,y):
    print 1 - Bio.Cluster.distancematrix((x,y), dist="k")[1][0]