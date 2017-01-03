'''
Created on Nov 15, 2016

@author: javi
'''

def mi(contact_map,x_nat_t,y_evol_t,x_nat_f,y_evol_f):
    import numpy as np
    import matplotlib.pyplot as plt
    N = 50
    #trues
    x_nat_t = [27.0,89.0,25.0]
    y_evol_t = [5.9,3.2,4.0]
    #falses
    x_nat_f = [3.0,1.0,3.0]
    y_evol_f = [1.1,0.32,1.0]
    
    
    
    #plt.scatter(x_nat_t, y_evol_t,color="b")
    #plt.scatter(x_nat_f, y_evol_f,color="r")
    x=[[1.0,2.4]]
    plt.scatter(x,color="r")
    
    plt.show()
    
mi()    