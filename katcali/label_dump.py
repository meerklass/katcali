
import numpy as np
from . import filter as kf
from . import diode as kd
from astropy.coordinates import SkyCoord
from astropy import units as u

def select_track(data,ant,pol):
    data.select(ants=ant,pol=pol,scans='track')
    scans_t=[]
    scans_tw=[]
    for s in data.scans():
        if data.shape[0]> 10: #50 modified @20250224
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
        if data.shape[0]> 10: #50 modified @20250224
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

def edge_dp_drop(dps,drop_num=1):
    dps=list(dps)
    dps.sort()
    dps=np.array(dps)
    l1=len(dps)
    dps=dps[drop_num:-drop_num]
    l2=len(dps)
    print ('# edge drop applied: '+str(l1)+' -> '+str(l2))
    return dps

def cal_dp_c1_split(dp_c1a,dp_c2a,dp_c3a,dp_c1b,dp_c2b,dp_c3b):
    i1=np.where(dp_c1a<np.min(dp_c2a))[0]
    dp_c1a1=dp_c1a[i1]
    i2=np.where(dp_c1a>np.max(dp_c3a))[0]
    dp_c1a3=dp_c1a[i2]
    dp_c1a2=dp_c1a[i1[-1]+1:i2[0]]
    assert(len(dp_c1a)==len(dp_c1a1)+len(dp_c1a2)+len(dp_c1a3))
    
    i1=np.where(dp_c1b<np.min(dp_c2b))[0]
    dp_c1b1=dp_c1b[i1]
    i2=np.where(dp_c1b>np.max(dp_c3b))[0]
    dp_c1b3=dp_c1b[i2]
    dp_c1b2=dp_c1b[i1[-1]+1:i2[0]]
    assert(len(dp_c1b)==len(dp_c1b1)+len(dp_c1b2)+len(dp_c1b3))
    
    
    dp_c1a1,dp_c1a2,dp_c1a3,dp_c1b1,dp_c1b2,dp_c1b3=list(dp_c1a1),list(dp_c1a2),list(dp_c1a3),list(dp_c1b1),list(dp_c1b2),list(dp_c1b3)
    dp_c1a1.sort()
    dp_c1a2.sort()
    dp_c1a3.sort()
    dp_c1b1.sort()
    dp_c1b2.sort()
    dp_c1b3.sort()

    return dp_c1a1,dp_c1a2,dp_c1a3,dp_c1b1,dp_c1b2,dp_c1b3

def cal_dp_c(fname,data,ant,pol,flags,ch,dp_tt,dp_ss,ang_deg, target_start=0, n_src_off=-1):
    sigma_level=10
    n_iter=3
    flags_1ch=flags[:,ch]
    dp_sb=dp_ss[0]
    dp_se=dp_ss[-1]
    data.select(ants=ant,pol=pol,scans='track',targets=target_start)
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

    data.select(ants=ant,pol=pol,scans='track',targets=target_start+1)
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

    if fname in ['1551055211','1551037708', '1579725085', '1580260015','1630519596','1631379874','1631387336','1631659886','1631667564'] or n_src_off==4:
        data.select(ants=ant,pol=pol,scans='track',targets=target_start+2)
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

        data.select(ants=ant,pol=pol,scans='track',targets=target_start+3)
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

        data.select(ants=ant,pol=pol,scans='track',targets=target_start+4)
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

        dp_c1a1,dp_c1a2,dp_c1a3,dp_c1b1,dp_c1b2,dp_c1b3=cal_dp_c1_split(dp_c1a,dp_c2a,dp_c3a,dp_c1b,dp_c2b,dp_c3b)
        dp_c1a1,dp_c1a2,dp_c1a3,dp_c1b1,dp_c1b2,dp_c1b3=edge_dp_drop(dp_c1a1),edge_dp_drop(dp_c1a2),edge_dp_drop(dp_c1a3),edge_dp_drop(dp_c1b1),edge_dp_drop(dp_c1b2),edge_dp_drop(dp_c1b3)
        dp_c1a=list(dp_c1a1)+list(dp_c1a2)+list(dp_c1a3)
        dp_c1b=list(dp_c1b1)+list(dp_c1b2)+list(dp_c1b3)
        
        dp_c0a, dp_c2a,dp_c3a,dp_c4a,dp_c0b,dp_c2b,dp_c3b,dp_c4b=edge_dp_drop(dp_c0a), edge_dp_drop(dp_c2a), edge_dp_drop(dp_c3a), edge_dp_drop(dp_c4a), edge_dp_drop(dp_c0b), edge_dp_drop(dp_c2b),edge_dp_drop(dp_c3b), edge_dp_drop(dp_c4b)
        #overwrite
        dp_ca=list(dp_c0a)+list(dp_c1a)+list(dp_c2a)+list(dp_c3a)+list(dp_c4a)
        dp_ca.sort()
        dp_cb=list(dp_c0b)+list(dp_c1b)+list(dp_c2b)+list(dp_c3b)+list(dp_c4b)
        dp_cb.sort()
    #result= dp_ca,dp_cb,dp_c0a, dp_c1a,dp_c2a,dp_c3a,dp_c4a,dp_c0b,dp_c1b,dp_c2b,dp_c3b,dp_c4b
    result= dp_ca,dp_cb,dp_c0a, dp_c1a,dp_c2a,dp_c3a,dp_c4a,dp_c0b,dp_c1b,dp_c2b,dp_c3b,dp_c4b,dp_c1a1,dp_c1a2,dp_c1a3,dp_c1b1,dp_c1b2,dp_c1b3
    data.select() #recover after select!!!
    data.select(ants=ant,pol=pol)
    return result

        
def cal_dp_u(dp_tt,dp_ss):
    dp_u=list(dp_tt)+list(dp_ss)
    dp_u.sort()
    dp_u=np.array(dp_u)
    return dp_u

##########################added for 2021 data#######################################

def cal_dp_label(data,flags,ant,pol,ch_ref,ang_deg):
    
    #### modified @20250224#####    
    dp_t,scans_tw=select_track(data,ant,pol)
    dp_s,scans_sw=select_scan(data,ant,pol)
    dp_w=select_waste(data,ant,pol) #includes slew,stop,scan_waste, track_waste     
    
    data.select() #reset
    data.select(ants=ant,pol=pol,scans='slew')
    dp_slew=data.dumps
    data.select(ants=ant,pol=pol,scans='stop')
    dp_stop=data.dumps
    ###data.select(ants=ant,pol=pol,scans='track')
    ###dp_t=data.dumps
    ###data.select(ants=ant,pol=pol,scans='scan')
    ###dp_s=data.dumps
    data.select() #reset
    data.select(ants=ant,pol=pol)

    ###dp_w=list(dp_slew)+list(dp_stop)
    ###dp_w.sort()
    
    #dp_sb=dp_s[0] #old version was dp_ss
    #dp_se=dp_s[-1]
        
    dp_tt=[]
    dp_ss=[]
    dp_f=[]
    flags_1ch=flags[:,ch_ref]
    for i in dp_t:
        if flags_1ch[i]==False and ang_deg[i]<=1:
            dp_tt.append(i)
        else:
            dp_f.append(i)

    for i in dp_s:
        if flags_1ch[i]==False:
            dp_ss.append(i)
        else:
            dp_f.append(i)

    return dp_tt,dp_ss,dp_f,dp_w, dp_t,dp_s,dp_slew,dp_stop
    #return dp_tt,dp_ss,dp_f,dp_w

    
def cal_label_intersec(dp_tt,dp_ss,nd_0,nd_1x):
    dp_tt,dp_ss,nd_0,nd_1x=np.array(dp_tt),np.array(dp_ss),np.array(nd_0),np.array(nd_1x)
    
    nd_t0=list(set(nd_0).intersection(set(dp_tt)))
    nd_t1x=list(set(nd_1x).intersection(set(dp_tt)))
    assert(len(nd_t0)+len(nd_t1x)==len(dp_tt))
    
    nd_s0=list(set(nd_0).intersection(set(dp_ss)))
    nd_s1x=list(set(nd_1x).intersection(set(dp_ss)))
    assert(len(nd_s0)+len(nd_s1x)==len(dp_ss))
    
    nd_t0,nd_t1x,nd_s0,nd_s1x=np.array(nd_t0),np.array(nd_t1x),np.array(nd_s0),np.array(nd_s1x)
    
    nd_t0_ca=nd_t0[nd_t0<dp_ss[0]]
    nd_t0_cb=nd_t0[nd_t0>dp_ss[-1]]
    assert(len(nd_t0)==len(nd_t0_ca)+len(nd_t0_cb))
    
    nd_t1x_ca=nd_t1x[nd_t1x<dp_ss[0]]
    nd_t1x_cb=nd_t1x[nd_t1x>dp_ss[-1]]
    assert(len(nd_t1x)==len(nd_t1x_ca)+len(nd_t1x_cb))
    
    nd_t0.sort(),nd_t1x.sort(),nd_s0.sort(),nd_s1x.sort(),nd_t0_ca.sort(),nd_t0_cb.sort(),nd_t1x_ca.sort(),nd_t1x_cb.sort()
    
    return nd_t0,nd_t1x,nd_s0,nd_s1x,nd_t0_ca,nd_t0_cb,nd_t1x_ca,nd_t1x_cb
    
    
def cal_label_intersec_complex(dp_tt,dp_ss,nd_0,nd_1x,nd_ratio):
    
    nd_t0,nd_t1x,nd_s0,nd_s1x,nd_t0_ca,nd_t0_cb,nd_t1x_ca,nd_t1x_cb=cal_label_intersec(dp_tt,dp_ss,nd_0,nd_1x) #load basic labels
    
    output=[]
    
    nd_ratio_group=kd.cal_nd_ratio_group(nd_ratio) #load ratio list
    print ('# nd ratio number in total: ' + str(len(nd_ratio_group))) 
    assert(nd_ratio_group[0]==0) # num0 is diode off so ignore it
    
    for i in range(1,len(nd_ratio_group)): #for each ratio number
        ratio_local=nd_ratio_group[i]
        print ('ratio'+str(i),end=" ")
        print (ratio_local,end="")
        print (': ',end="")
        
        label0=np.where(abs(nd_ratio-ratio_local)<1e-5)[0] #all dump with this ratio
        label0=list(set(label0).intersection(set(dp_ss+dp_tt))) #delete the waste dumps
        
        #for RFI flagging we need distinguish track-I, scan, track-II
        label1=list(set(label0).intersection(set(nd_t1x_ca)))
        label2=list(set(label0).intersection(set(nd_s1x)))
        label3=list(set(label0).intersection(set(nd_t1x_cb)))
        label1.sort(),label2.sort(),label3.sort()
        print (len(label0),len(label1),len(label2),len(label3))
        assert(len(label1)+len(label2)+len(label3)==len(label0))
        output.append(label1)
        output.append(label2)
        output.append(label3)
    #l=np.shape(output)[0]
    l=len(output)
    assert(l==3*(len(nd_ratio_group)-1))
    print ('# total label groups for nd_1x: '+str(l))
    return output
    
def cal_ptr_mask(p,p_radec,nd_s0, dp_sb,dp_se,ang_lim):
    dp_ptr_list=[]

    for i in range(len(p_radec)):
        #print i
        p_ra,p_dec=p_radec[i]
        c = SkyCoord(p_ra*u.deg,  p_dec*u.deg, frame='icrs')
        #print c 
        p_ang=(c.separation(p)/u.deg)[:,0]
        #print p_ang
        dp_l=np.where(p_ang<ang_lim)[0]
        #print dp_l
        for j in range(len(dp_l)):
            if dp_l[j]>dp_sb and dp_l[j]<=dp_se:
                dp_ptr_list.append(dp_l[j])

    list(set(dp_ptr_list))
    dp_ptr_list.sort()
    #dp_s0_ptr=list(set(dp_ptr_list).intersection(nd_s0))
    #dp_s0_ptr.sort()
    return dp_ptr_list
