import numpy as np
from scipy.interpolate import UnivariateSpline
from astropy.stats import sigma_clip

def deg_filter(dump, ang_deg_list, sigma_level, n_iter):
    ang_mean=ang_deg_list[dump].mean()
    ang_std=ang_deg_list[dump].std()
    dump_iter=list(dump)
    #print 'deg_filter'
    print 'deg filter start: '+str(ang_mean)+'+/-'+str(ang_std)
    
    for n in range(n_iter):
        for i in dump_iter:
            bias=abs(ang_deg_list[i]-ang_mean)
            if bias > sigma_level*ang_std:
                #print i
                dump_iter.remove(i)
        ang_mean=ang_deg_list[dump_iter].mean()
        ang_std=ang_deg_list[dump_iter].std()
        #print str(n)+': '+str(ang_mean)+'+/-'+str(ang_std)
    print 'deg filter end: '+str(ang_mean)+'+/-'+str(ang_std)+'\n'
    return dump_iter

def curve_filter(x,y,sigma=4,k=5):
    a=[]
    b=[]
    for i in range(len(x)):
        if y[i] is not None and np.isnan(y[i]) is not True: 
            a.append(x[i])
            b.append(y[i])
    a=np.array(a)
    b=np.array(b)
    l=UnivariateSpline(a,b,k=k,s=None)(a)
    filtered_data = sigma_clip(b-l, sigma=sigma, iters=5)
    return a[filtered_data.mask==False] ###output is what want to keep
