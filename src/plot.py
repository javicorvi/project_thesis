'''
Created on Jan 12, 2017

@author: javi
'''
import itertools
import matplotlib.pyplot as plt
import numpy as np
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
    plt.show()
    plt.gcf().clear()
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
    plt.show()   
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
    plt.show()
    plt.gcf().clear()
'''
Test function not in use
'''
def contact_map_(contact_map, outputpath):
    np.set_printoptions(threshold=np.nan)
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
    x_nat_t = [27.0,12.0,25.0]
    y_nat_t = [35.0,101.0,37]
    
    x_evol_t = [35.0,101.0,37]
    y_evol_t = [27.0,12.0,25.0]
    plt.scatter(x_nat_t, y_nat_t,color="b")
    plt.scatter(x_evol_t, y_evol_t,color="r")
    g = np.floor(cmap)
    print len(cmap)
    plt.imshow(g, cmap='Greys')   
    #plt.show()
    plt.savefig(outputpath + "/contact_map_mi.png")
    plt.gcf().clear()
