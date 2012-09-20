import numpy

def errfunc(P, x, y, model):                                                   
    return y-model(P, x) 

def gaussian(G=[0,1,1], x=1):                                                  
    return G[2]*numpy.exp(-0.5*((x-G[0])/G[1])**2) 
