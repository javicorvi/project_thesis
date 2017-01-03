import Bio.Cluster
import utils  as utils
import os
import time
from plot import auc as plot
from sklearn import preprocessing
def spearman_zmip(zmip_natural_result_path, mi_results_path, pattern_array,contact_map_path,outputpath):
    #levanto el zmip natural
    zmip_natural = utils.load_zmip(zmip_natural_result_path)
    utils.order(zmip_natural)
    for filename in os.listdir(mi_results_path):
        if filename.endswith(".dat") & any(r in filename for r in pattern_array):
            
            print " Calculation of : " + filename + " with contact map " + contact_map_path
            #levanto el zmip evolucionado 
            zmip_evol = utils.load_zmip(mi_results_path + filename)
            
            #sincronizo las senales de coevolucion para calcular spearman rank correlation
            m,m2=utils.sincronice_mi(zmip_natural, zmip_evol)
            m_=[row [2] for row in m]
            m2_=[row[2] for row in m2]
            value_spearman = spearman(m_,m2_)
            print "Spearman total " + str(value_spearman)
            
            #MI PLOT
            contact_map=utils.load_contact_map(contact_map_path)
            
            plot.contact_map(contact_map)
            #Test load_contact_map
            #v=contact_map[4][2]
            #v=contact_map[4][1]
            #v=contact_map[4][3]
            #v=contact_map[12][4]
            import numpy as np        
            m1_norm,m2_norm=normalice_(m_,m2_)
            m_np = np.c_[ np.asarray(m), np.asarray(m1_norm) ]
            m2_np = np.c_[ np.asarray(m2), np.asarray(m2_norm) ]   
            x_nat_t = []
            y_evol_t = []
            x_nat_f = []
            y_evol_f = [] 
            for x, y in map(None, m_np, m2_np):
                #se le resta a uno porque los arrays comienzan en la posicion 0
                pos1 = int(x[0]-1)
                pos2 = int(x[1]-1)
                v = contact_map[pos1][pos2]
                if( 'true' in v ):
                    x_nat_t.append(x[3])
                    y_evol_t.append(y[3])
                else:
                    x_nat_f.append(x[3])
                    y_evol_f.append(y[3])
                
            
            
            plot.mi(x_nat_t,y_evol_t,x_nat_f,y_evol_f,outputpath+filename+'.png',filename)
            
            #ordeno zmip evolucionado sincronizado 
            utils.order(m2)
            
            print "TOTAL PAR POSITIONS " + str(len(zmip_natural))
            top_rank(zmip_natural,m2,1)
            top_rank(zmip_natural,m2,5)
            top_rank(zmip_natural,m2,10)
            top_rank(zmip_natural,m2,12)
            print '************************************************************************'
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
                j.append(e[2])   
                data.append(j)
                break    
    return  data
def top_rank(x,y,top):
    num = len(x)*top/100
    a=x[0:num]
    b=y[0:num]
    data=matches_coevolved_positions(a,b)
    print "TOP "  + str(top) + "% PAR POSITIONS " + str(num)
    print "HITS " + str(len(data))
    #print data
def kendall(x,y):
    z = [1.65, 2.64, 2.64, 6.95]
    print 1 - Bio.Cluster.distancematrix((x,y), dist="k")[1][0]

    
        
    