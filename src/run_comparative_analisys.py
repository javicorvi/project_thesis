__author__ = 'javicorvi'
# coding=utf-8
import dataanalisys
import util
import numpy as np 
import plot
data_paths=[
            "../2trx_s3_w3/mi_data/zmipsequences_2trx_edit-beta1.00-nsus3.00-runs20000.fasta_0.62.cluster.dat",
            "../1thx_s3_w3/mi_data/zmipsequences_1thx_edit-beta1.00-nsus3.00-runs20000.fasta_0.62.cluster.dat"]

contact_map_path="../joined_s3w0/contact_map/contact_map_2trx_1thx.dat"
filename="joined"
result_file="joined.dat"
mi1=util.load_zmip(data_paths[0])
mi2=util.load_zmip(data_paths[1])
minatural=util.load_zmip("../2trx_s3_w3/natural/zmip_PF00085_THIO_ECOLI_reference_1thx.dat")



m,m1=util.sincronice_mi(minatural, mi1)
m,m2=util.sincronice_mi(minatural, mi2)
m_=[row [2] for row in m]
m1_=[row[2] for row in m1]
m2_=[row[2] for row in m2]

m_norm,m1_norm,m2_norm=dataanalisys.normalice_desarrollo(m_,m1_,m2_)
m_np = np.c_[ np.asarray(m), np.asarray(m_norm) ]
m1_np = np.c_[ np.asarray(m1), np.asarray(m1_norm) ]   
m2_np = np.c_[ np.asarray(m2), np.asarray(m2_norm) ]  
#array for contact mi matrix comparission
x_nat_t = []
y_evol_t = []
y_evol2_t = []
x_nat_f = []
y_evol_f = []
y_evol2_f = []

contact_map=util.load_contact_map(contact_map_path)

#scores mi, for roc_curve 
scores_nat = []
scores_evol = []
scores_evol2 = []
#1 contact
#0 no contact
y_true = []
for x, y, f in map(None, m_np, m1_np,m2_np):
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
    scores_evol2.append(f[3])
    if(v == 0):
        x_nat_f.append(x[3])
        y_evol_f.append(y[3])
        y_evol2_f.append(f[3])
        y_true.append(0)
    else:
        x_nat_t.append(x[3])
        y_evol_t.append(y[3])
        y_evol2_t.append(f[3])
        y_true.append(1)
            
labels=['Natural', 'Evol 2TRX','Evol 1THX']
scores=[]
scores.append(scores_nat)
scores.append(scores_evol)
scores.append(scores_evol2)
#plot.roc_curve(y_true,scores_nat,scores_evol)
colors = ['blue', 'red', 'green']
plot.roc_curve(y_true,scores,labels,colors, filename+'_roc_curve.png')
plot.contacts_with_mi_desarrollo(x_nat_t,y_evol_t,y_evol2_t,x_nat_f,y_evol_f,y_evol2_f,filename+'contacts_with_mi.png',filename)

 #dataanalisys.top_rank_desa(minatural,mi,mi2,1,contact_map,"plot_joined_miÇ_top_1percent_withcon.png",filename,result_file)
