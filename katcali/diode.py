import numpy as np
from scipy import stats

def label_nd_injection(fname,vis, timestamps, dp_ss, dump_period):

    ######set a jump limit###############
    f=10.
    if fname in ['1555793534','1551055211','1551037708','1555775533']:
        f=10
    if fname=='1556120503':
        f=2
    if fname=='1556052116':
        f=15.

    ch_plot0=800 #only for edge detection


    lmax=abs(np.nanmax(vis[dp_ss,ch_plot0]))
    lmin=abs(np.nanmin(vis[dp_ss,ch_plot0]))
    lim=(lmax-lmin)/f

    #print lmax,lmin,lim

    mark=[]
    nd_1=[]
    nd_1a=[]
    nd_1b=[]
    for i in range(1,len(timestamps)):
        if (np.abs(vis[i,ch_plot0])-np.abs(vis[i-1,ch_plot0]) >lim #have a jump
            and vis[i,ch_plot0]>0 and vis[i-1,ch_plot0]> 0
            and timestamps[i-1]-timestamps[0] -dump_period/2. not in mark): #not jump the one before
        
            m=timestamps[i]-timestamps[0]-dump_period/2.
            #print i,m
            mark.append(m) #for plot
            nd_1.append(i) #on
            nd_1a.append(i) #on1
            if i+1<len(timestamps):
                nd_1.append(i+1) #on
                nd_1b.append(i+1) #on2
    return mark,nd_1,nd_1a,nd_1b,lmin,lmax

                
def gap_list(list):
    gap_list=[]
    for i in range(1,len(list)):
        gap_list.append(list[i]-list[i-1])
    
    gap_list=np.array(gap_list)
    return gap_list

def gap_mode(list):
    gap_list=[]
    for i in range(1,len(list)):
        gap_list.append(list[i]-list[i-1])
    
    gap_list=np.array(gap_list)
    mode=stats.mode(gap_list)[0][0]
    return mode

def cal_nd_wro_0(nd_1a):
    nd_gap_list=gap_list(nd_1a)
    nd_gap_mode=gap_mode(nd_1a)
    nd_wro_0=np.where(nd_gap_list != nd_gap_mode)[0]
    return nd_gap_mode, nd_wro_0

def cal_t_line(fname, timestamps,nd_set, nd_cycle, dump_period):
    t_line=[]
    nd_tb=nd_set-timestamps[0]
    nd_te=dump_period*len(timestamps)

    for t in np.arange(nd_tb, nd_te+dump_period, nd_cycle ) :
        t_line.append(t)
    #######below are caused by the diode lead time###############
    if fname=='1551037708':
        t_line=t_line[5:]
    
    if fname=='1551055211':
        t_line=t_line[6:]
    
    if fname=='1553966342':
        t_line=t_line[2:]
    
    if fname=='1556034219':
        t_line=t_line[2:]
    
    if fname=='1554156377':
        t_line=t_line[1:]
    
    if fname=='1556138397':
        t_line=t_line[2:-1]
        
    if fname=='1556052116':
        t_line=t_line[2:]
    
    if fname=='1555775533':
        t_line=t_line[2:]
    
    if fname=='1556120503':
        t_line=t_line[2:]
    
    if fname=='1555793534':
        t_line=t_line[2:]
    
    if fname=='1561650779':
        t_line=t_line[3:-1]

    if fname=='1562857793':
        t_line=t_line[3:]
    ########################################################################
    return t_line
    
def call_nd_1a_param(fname):
    nd_1a_gap=10
    nd_1a0=-999
    if fname=='1551055211':
        nd_1a0=1
    if fname=='1555793534':
        nd_1a0 =8
    if fname=='1551037708':
        nd_1a0 =4
    if fname=='1553966342':
        nd_1a0 =4
    if fname=='1554156377':
        nd_1a0 =5
    if fname=='1555775533':
        nd_1a0 =0
    if fname=='1556034219':
        nd_1a0 =-1 #to keep nd_1b in the list
    if fname=='1556120503':
        nd_1a0 =0
    if fname=='1556138397':
        nd_1a0 =2
    if fname=='1556052116':
        nd_1a0 =1
    if fname=='1561650779':
        nd_1a0 =7
    if fname=='1562857793':
        nd_1a0 =7    
    if nd_1a0==-999:
        print('no record, can ask astro.jywang@gmail.com')
    return nd_1a_gap,nd_1a0

def call_nd_1_list(fname,timestamps):
    nd_1a_gap,nd_1a0=call_nd_1a_param(fname)
    nd_1aa=[]
    for i in range(1000):
        a=nd_1a0+i*nd_1a_gap
        if a >= 0:
            if a >= len(timestamps):
                break
            nd_1aa.append(a)
          
    nd_1bb=[]
    for i in range(1000):
        a=nd_1a0+1+i*nd_1a_gap
        if a>=0:
            if a >= len(timestamps):
                break
            max=i
            nd_1bb.append(a)
    nd_11=list(nd_1aa)+list(nd_1bb)
    nd_11.sort()
    
    nd_00=[]
    for i in range(len(timestamps)):
        if i not in nd_11:
            nd_00.append(i)

    assert(len(nd_00)+len(nd_11)==len(timestamps))
    return nd_1aa,nd_1bb,nd_11,nd_00
       
def cal_nds_list(dp_ss,nd_1a,nd_1b,nd_1,nd_0): #dp_ss/dp_s
    nd_s0=[]
    for i in dp_ss:
        if i in nd_0:
            nd_s0.append(i)
    #print np.shape(nd_s0)       

    nd_s1=[]
    nd_s1a=[]
    nd_s1b=[]
    for i in dp_ss:
        if i in nd_1:
            nd_s1.append(i)
        if i in nd_1a:
            nd_s1a.append(i)
        if i in nd_1b:
            nd_s1b.append(i)
    return nd_s1a,nd_s1b,nd_s1,nd_s0

def cal_ndt_list(dp_tt,nd_1a,nd_1b,nd_1,nd_0):#dp_tt/dp_s
    nd_t0=[]
    for i in dp_tt:
        if i in nd_0:
            nd_t0.append(i)
    #print np.shape(nd_t0)       

    nd_t1=[]
    nd_t1a=[]
    nd_t1b=[]
    for i in dp_tt:
        if i in nd_1:
            nd_t1.append(i)
        if i in nd_1a:
            nd_t1a.append(i)
        if i in nd_1b:
            nd_t1b.append(i)
    return nd_t1a,nd_t1b,nd_t1,nd_t0

def cal_az_edge_list(az,dp_sb,dp_se):
    dp_az_min=[]
    dp_az_max=[]
    for i in range(dp_sb,dp_se+1): #dp_w included
        if az[i]<az[i-1] and az[i]<az[i+1]:
            dp_az_min.append(i)
        if az[i]>az[i-1] and az[i]>az[i+1]:  
            dp_az_max.append(i)
    dp_az_edge=list(dp_az_min+dp_az_max)
    dp_az_edge.sort()
    return dp_az_min,dp_az_max,dp_az_edge
