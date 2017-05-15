'''
Created on Jan 12, 2017

@author: javi
'''
import itertools
import matplotlib.pyplot as plt
import numpy as np
from sklearn import metrics
'''
Generates a plot describe the contacts with de MI of the natural and evolution msa
'''
def contacts_with_mi(x_nat_t,y_evol_t,x_nat_f,y_evol_f,output_path,filename):
    #trues
    #x_nat_t = [27.0,89.0,25.0]
    #y_evol_t = [5.9,3.2,4.0]
    #falses
    #x_nat_f = [3.0,1.0,3.0]
    #y_evol_f = [1.1,0.32,1.0]
    #plt.scatter(x_nat_t, y_evol_t,color="b")
    #plt.scatter(x_nat_f, y_evol_f,color="r")
    plt.axis([0, 1, 0, 1])
    no = plt.scatter(y_evol_f, x_nat_f,color="r")
    co = plt.scatter(y_evol_t, x_nat_t,color="b")
    plt.legend((no, co),
           ('No Contact', 'Contact'),
           scatterpoints=1,
           loc='lower left',
           ncol=3,
           fontsize=8)
    #plt.set_title(filename)
    plt.ylabel('Natural MI Values')
    plt.xlabel('Evolution MI Values')
    plt.title(filename)
    
    #x=[[1.0,2.4]]
    #plt.scatter(x,color="r")
    plt.savefig(output_path)
    #plt.show()
    plt.gcf().clear()

'''
Deprecated
'''
def roc_curve_dep(y_true,scores_nat,scores_evol):
    fpr, tpr, _ = metrics.roc_curve(y_true, scores_nat)
    roc_auc = metrics.auc(fpr, tpr)
    
    fpr2, tpr2, _ = metrics.roc_curve(y_true, scores_evol)
    roc_auc2 = metrics.auc(fpr2, tpr2)      
            
    plt.figure()
    lw = 2
    
    
    plt.plot(fpr, tpr, color='blue', lw=lw, label='ROC curve Nat (area = %0.2f)' % roc_auc)
    plt.plot(fpr2, tpr2, color='red', lw=lw, label='ROC curve Evol (area = %0.2f)' % roc_auc2)
    
    plt.plot([0, 1], [0, 1], color='black', lw=lw, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver operating characteristic example')
    plt.legend(loc="lower right")
    plt.show()
    plt.gcf().clear()
    
def roc_curve(df,index,y_true,scores,labels,colors,output_file):
    fpr = dict()
    tpr = dict()
    roc_auc = dict()
    lw = 2
    plt.figure()
    #for i, (a, b) in enumerate(zip(alist, blist)):
    for i,(label,color) in enumerate(zip(labels, colors)):
        fpr[i], tpr[i], _ = metrics.roc_curve(y_true, scores[i])
        partial_auc_value_0_1 = partial_auc(fpr[i], tpr[i], 0.1)
        auc = metrics.auc(fpr[i], tpr[i])
        roc_auc[i] = auc
        if(i==0):
            df.set_value(index,'auc_nat',auc)
            df.set_value(index,'auc_nat_01',partial_auc_value_0_1)
        else:
            df.set_value(index,'auc',auc)
            df.set_value(index,'auc_01',partial_auc_value_0_1)
        plt.plot(fpr[i], tpr[i], color=color, lw=lw,label='ROC curve of class {0} (auc = {1:0.2f} | auc 0.1 = {2:0.2f})'''.format(label, roc_auc[i], partial_auc_value_0_1))
    plt.plot([0, 1], [0, 1], color='black', lw=lw, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC CONTACTS MI')
    plt.legend(loc="lower right")
    plt.savefig(output_file)
    plt.show()
    plt.gcf().clear()

def partial_auc(fpr, tpr, max_fpr):
    idx = np.where(fpr <= max_fpr)[0]
    # linearly interpolate the ROC curve until max_fpr
    idx_last = idx.max()
    idx_next = idx_last + 1
    xc = [fpr[idx_last], fpr[idx_next]]
    yc = [tpr[idx_last], fpr[idx_next]]
    tpr = np.r_[tpr[idx], np.interp(max_fpr, xc, yc)]
    fpr = np.r_[fpr[idx], max_fpr]
    partial_roc = metrics.auc(fpr, tpr, reorder=True)
    # standardize result to lie between 0.5 and 1
    min_area = max_fpr**2/2
    max_area = max_fpr
    return 0.5*(1+(partial_roc-min_area)/(max_area-min_area))    
'''
Generate the plot for the diferents auc taking into account the beta, nsus and runs
Open the result_auc_path parse all the results and plot them into 4 subplots
one for each runs
'''
def plot_auc(result_auc_path,output_path,beta,runs,nsus):
    ys_run = np.zeros((len(runs),len(beta), len(nsus)))
    for id_r,r in enumerate(runs):  
        for id_b,b in enumerate(beta):   
            #data_beta = [] * len(nsus)
            for id_n,n in enumerate(nsus):    
                source = open(result_auc_path, 'r')
                for line in source:
                    if (('beta'+b in line) & ('nsus'+n in line) & ('runs'+r+'.' in line)):
                        str,value=line.split(",")
                        #data_beta[id_n]=value
                        ys_run[id_r,id_b,id_n]=value
                source.close()  
    print ys_run  
    
    f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2 , sharex='col', sharey='row')
    #f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2)
    colors = itertools.cycle([ "yellow","orange" ,"green", "red", "blue",  "black"])
    subplt = [ax1,ax2,ax3,ax4]
    
    axis=[0, 6, 0, 1.0]
    ax1.axis(axis)
    ax2.axis(axis)
    ax3.axis(axis)
    ax4.axis(axis)
    
    #for splt in subplt:
    for id_splt,splt in enumerate(subplt):
        i=0
        splt.set_title('Runs ' + runs[id_splt])
        
        for y, c in zip(ys_run[id_splt], colors):
            splt.scatter(nsus, y, color=c,picker=True)  
            splt.plot(nsus, y,color=c, label=beta[i])
            i=i+1    
    
    
    # Now add the legend with some customizations.
    
    f.text(0.5, 0.04, 'NSUS', ha='center')
    f.text(0.04, 0.5, 'AUC', va='center', rotation='vertical')
    #f.canvas.mpl_connect('pick_event', onpick3)
    plt.legend(loc=1,prop={'size':10},title="Betas")
    #plt.show()   
    #plt.savefig(output_path+'/auc.png')
    plt.gcf().clear()

'''

'''
def contact_map_with_top_rank_mi(contact_map, x_nat, y_nat, x_evol,y_evol,output_path,filename):
    #np.set_printoptions(threshold=np.nan)
    #print contact_map
    #contact_map[contact_map == 'false']=0
    #contact_map[contact_map == 'true']=1
    #cmap=np.array(contact_map, dtype='i4')
    #print cmap2
    #print len(cmap2)
    #cmap=np.random.randint(2, size=(108, 108))
    #print len(cmap)
    #np.random.seed(101)
    cmap=contact_map
    #x_nat_t = [27.0,89.0,25.0]
    #y_evol_t = [5.9,3.2,4.0]
    plt.scatter(x_nat, y_nat,color="b",s=40,  marker=(5, 2))
    plt.scatter(y_evol, x_evol, color="r",s=40,  marker=(5, 2))
    g = np.floor(cmap)
    plt.imshow(g, cmap='Greys')
    plt.title(filename)
    plt.savefig(output_path)
    #plt.show()
    plt.gcf().clear()

def contact_map_with_top_rank_mi_desarrollo(contact_map, x_nat, y_nat, x_evol,y_evol,x_evol2,y_evol2,output_path,filename):
    cmap=contact_map
    plt.scatter(x_nat, y_nat,color="b",s=40,  marker=(5, 2))
    plt.scatter(y_evol, x_evol, color="r",s=40,  marker=(5, 2))
    plt.scatter(y_evol2, x_evol2, color="g",s=40,  marker=(5, 2))
    g = np.floor(cmap)
    plt.imshow(g, cmap='Greys')
    plt.title(filename)
    plt.savefig(output_path)
    #plt.show()
    plt.gcf().clear()

'''
Plot the contact map and saved in output_file
'''
def contact_map(contact_map, output_file):
    import matplotlib.cm as cm
    '''
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    from matplotlib.colors import LogNorm
    import numpy as np
    x, y, z = np.loadtxt('data.txt', unpack=True)
    '''
    cmap=contact_map
    
    plt.imshow(cmap,cmap=cm.Blues)
    plt.colorbar()
    #plt.imshow(cmap,cmap=cm.hot)
    plt.savefig(output_file)
    plt.show()
    plt.gcf().clear()
'''
Generates a plot describe the contacts with de MI of the natural and evolution msa
'''
def contacts_with_mi_desarrollo(x_nat_t,y_evol_t,y_evol2_t,x_nat_f,y_evol_f,y_evol2_f,output_path,filename):
    #trues
    #x_nat_t = [27.0,89.0,25.0]
    #y_evol_t = [5.9,3.2,4.0]
    #falses
    #x_nat_f = [3.0,1.0,3.0]
    #y_evol_f = [1.1,0.32,1.0]
    #plt.scatter(x_nat_t, y_evol_t,color="b")
    #plt.scatter(x_nat_f, y_evol_f,color="r")
    plt.axis([0, 1, 0, 1])
    no = plt.scatter(y_evol_f, x_nat_f,color="r")
    no2= plt.scatter(y_evol2_f, x_nat_f,color="r")
    co = plt.scatter(y_evol_t, x_nat_t,color="b")
    co2 = plt.scatter(y_evol2_t, x_nat_t,color="g")
    plt.legend((no, co,co2),
           ('No Contact', 'Contact 2TRX', 'CONTACT 1THX'),
           scatterpoints=1,
           loc='lower left',
           ncol=3,
           fontsize=8)
    #plt.set_title(filename)
    plt.ylabel('Natural MI Values')
    plt.xlabel('Evolution MI Values')
    plt.title(filename)
    
    #x=[[1.0,2.4]]
    #plt.scatter(x,color="r")
    plt.savefig(output_path)
    plt.show()
    plt.gcf().clear()
'''
Plot the conservation of the MSA evolutionated and the Natural MSA
'''
def conservation_between_msas(msas_entropy, output_file,natural_line_style='-'):
    for index,msa_entropy in enumerate(msas_entropy):
        if(index==0):
            plt.plot(msa_entropy[0],color=np.random.rand(3,1), label=msa_entropy[1],linewidth=5,linestyle=natural_line_style)
        else:
            plt.plot(msa_entropy[0],color=np.random.rand(3,1), label=msa_entropy[1])
    plt.ylabel('Bits Entropy')
    plt.xlabel('Position')
    plt.legend(loc=1,prop={'size':10},title="PDB")
    plt.savefig(output_file)
    plt.show()  
    plt.gcf().clear()

