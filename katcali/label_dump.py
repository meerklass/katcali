
import numpy as np
from . import filter as kf

def select_track(data,ant,pol):
    data.select(ants=ant,pol=pol,scans='track')
    scans_t=[]
    scans_tw=[]
    for s in data.scans():
        if data.shape[0]> 50:
            scans_t.append(data.scan_indices[0])
        else:
            scans_tw.append(data.scan_indices[0])
    data.select(ants=ant,pol=pol,scans=scans_t)
    dp_t=data.dumps
    data.select() #recover after select!!!
    data.select(ants=ant,pol=pol)
    return dp_t,scans_tw

def select_scan(data,ant,pol):
    data.select(ants=ant,pol=pol,scans='scan')
    scans_s=[]
    scans_sw=[]
    for s in data.scans():
        if data.shape[0]> 50:
            scans_s.append(data.scan_indices[0])
        else:
            scans_sw.append(data.scan_indices[0])
    data.select(ants=ant,pol=pol,scans=scans_s)
    dp_s=data.dumps
    #dp_sb=data.dumps[0]
    #dp_se=data.dumps[-1]
    data.select() #recover after select!!!
    data.select(ants=ant,pol=pol)
    return dp_s,scans_sw

def select_waste(data,ant,pol):
    data.select(ants=ant,pol=pol,scans=('slew','stop'))
    dp_w1=data.dumps
    
    dp_t,scans_tw=select_track(data,ant,pol)
    dp_s,scans_sw=select_scan(data,ant,pol)
    
    data.select(ants=ant,pol=pol,scans=scans_tw+scans_sw)
    dp_w2=data.dumps
    dp_w=list(dp_w1)+list(dp_w2)
    dp_w.sort()
    dp_w=np.array(dp_w)
    data.select() #recover after select!!!
    data.select(ants=ant,pol=pol)
    return dp_w

def label_dump_1ch(data,ant,pol,flags,ch):
    dp_t,scans_tw=select_track(data,ant,pol)
    dp_s,scans_sw=select_scan(data,ant,pol)
    dp_w=select_waste(data,ant,pol)
    
    flags_1ch=flags[:,ch]

    dp_tt=[]
    dp_ss=[]
    dp_f=[]

    for i in dp_t:
        if flags_1ch[i]==False:
            dp_tt.append(i)
        else:
            dp_f.append(i)
        
    for i in dp_s:
        if flags_1ch[i]==False:
            dp_ss.append(i)
        else:
            dp_f.append(i)
    return dp_tt,dp_ss,dp_f,dp_t,dp_s
    
def cal_dp_c(fname,data,ant,pol,flags,ch,dp_tt,dp_ss,ang_deg):
    sigma_level=10
    n_iter=3
    flags_1ch=flags[:,ch]
    dp_sb=dp_ss[0]
    dp_se=dp_ss[-1]
    data.select(ants=ant,pol=pol,scans='track',targets=0)
    dp_c0=data.dumps
    for i in dp_c0:
        if i not in dp_tt or flags_1ch[i]==True:
            dp_c0=list(dp_c0)
            dp_c0.remove(i)
            #print 'rm '+str(i)+' from dp_c0'
    dp_c0=kf.deg_filter(dp_c0,ang_deg,sigma_level,n_iter)
    dp_c0=np.array(dp_c0)
    dp_c0a=dp_c0[dp_c0<dp_sb]
    dp_c0b=dp_c0[dp_c0>dp_se]

    data.select(ants=ant,pol=pol,scans='track',targets=1)
    dp_c1=data.dumps
    for i in dp_c1:
        if i not in dp_tt or flags_1ch[i]==True:
            dp_c1=list(dp_c1)
            dp_c1.remove(i)
            #print 'rm '+str(i)+' from dp_c1'
    dp_c1=kf.deg_filter(dp_c1,ang_deg,sigma_level,n_iter)
    dp_c1=np.array(dp_c1)
    dp_c1a=dp_c1[dp_c1<dp_sb]
    dp_c1b=dp_c1[dp_c1>dp_se]

    dp_ca=list(dp_c0a)+list(dp_c1a)#+list(dp_c2a)+list(dp_c3a)+list(dp_c4a)
    dp_ca.sort()
    dp_cb=list(dp_c0b)+list(dp_c1b)#+list(dp_c2b)+list(dp_c3b)+list(dp_c4b)
    dp_cb.sort()
    result= dp_ca,dp_cb,dp_c0a, dp_c1a,dp_c0b,dp_c1b
    
    #######################all have above#######################

    if fname in ['1551055211','1551037708', '1579725085', '1580260015']:
        data.select(ants=ant,pol=pol,scans='track',targets=2)
        dp_c2=data.dumps
        for i in dp_c2:
            if i not in dp_tt or flags_1ch[i]==True:
                dp_c2=list(dp_c2)
                dp_c2.remove(i)
                #print 'rm '+str(i)+' from dp_c2'
        dp_c2=kf.deg_filter(dp_c2,ang_deg,sigma_level,n_iter)
        dp_c2=np.array(dp_c2)
        dp_c2a=dp_c2[dp_c2<dp_sb]
        dp_c2b=dp_c2[dp_c2>dp_se]

        data.select(ants=ant,pol=pol,scans='track',targets=3)
        dp_c3=data.dumps
        for i in dp_c3:
            if i not in dp_tt or flags_1ch[i]==True:
                dp_c3=list(dp_c3)
                dp_c3.remove(i)
                #print 'rm '+str(i)+' from dp_c3'
        dp_c3=kf.deg_filter(dp_c3,ang_deg,sigma_level,n_iter)
        dp_c3=np.array(dp_c3)
        dp_c3a=dp_c3[dp_c3<dp_sb]
        dp_c3b=dp_c3[dp_c3>dp_se]

        data.select(ants=ant,pol=pol,scans='track',targets=4)
        dp_c4=data.dumps
        for i in dp_c4:
            if i not in dp_tt or flags_1ch[i]==True:
                dp_c4=list(dp_c4)
                dp_c4.remove(i)
                #print 'rm '+str(i)+' from dp_c4'
        dp_c4=kf.deg_filter(dp_c4,ang_deg,sigma_level,n_iter)
        dp_c4=np.array(dp_c4)
        dp_c4a=dp_c4[dp_c4<dp_sb]
        dp_c4b=dp_c4[dp_c4>dp_se]
        #overwrite
        dp_ca=list(dp_c0a)+list(dp_c1a)+list(dp_c2a)+list(dp_c3a)+list(dp_c4a)
        dp_ca.sort()
        dp_cb=list(dp_c0b)+list(dp_c1b)+list(dp_c2b)+list(dp_c3b)+list(dp_c4b)
        dp_cb.sort()
        result= dp_ca,dp_cb,dp_c0a, dp_c1a,dp_c2a,dp_c3a,dp_c4a,dp_c0b,dp_c1b,dp_c2b,dp_c3b,dp_c4b
    data.select() #recover after select!!!
    data.select(ants=ant,pol=pol)
    return result

        
def cal_dp_u(dp_tt,dp_ss):
    dp_u=list(dp_tt)+list(dp_ss)
    dp_u.sort()
    dp_u=np.array(dp_u)
    return dp_u
