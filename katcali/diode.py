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
    if fname=='1579725085':
        nd_1a0 =3
    if fname=='1580260015':
        nd_1a0 =2
        
    if nd_1a0==-999:
        print ('no record, can ask astro.jywang@gmail.com')
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

#label track dumps as before and after scan
def cal_nd_t_c_list(nd_t0,nd_t1a, nd_t1b, dp_sb,dp_se):
    
    nd_t0_ca,nd_t0_cb,nd_t1a_ca,nd_t1a_cb,nd_t1b_ca,nd_t1b_cb=[],[],[],[],[],[]
    
    for i in nd_t0:
        if i< dp_sb:
            nd_t0_ca.append(i)
        if i> dp_se:
            nd_t0_cb.append(i)
    for i in nd_t1a:
        if i< dp_sb:
            nd_t1a_ca.append(i)
        if i> dp_se:
            nd_t1a_cb.append(i)
    for i in nd_t1b:
        if i< dp_sb:
            nd_t1b_ca.append(i)
        if i> dp_se:
            nd_t1b_cb.append(i)
    nd_t0_ca,nd_t0_cb,nd_t1a_ca,nd_t1a_cb,nd_t1b_ca,nd_t1b_cb=np.array(nd_t0_ca),np.array(nd_t0_cb),np.array(nd_t1a_ca),np.array(nd_t1a_cb),np.array(nd_t1b_ca),np.array(nd_t1b_cb)
    return nd_t0_ca,nd_t0_cb,nd_t1a_ca,nd_t1a_cb,nd_t1b_ca,nd_t1b_cb


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

##########################added for 2021 data#######################################

def cal_nd_basic_para(fname):
    
    nd_on_time=-999
    ###UHF#####

    if fname=='1684087370':
        nd_on_time=0.584099237647
        nd_cycle=19.4699745882
        nd_set=1684087370.34

    if fname=='1675632179':
        nd_on_time=0.584099237647
        nd_cycle=19.4699745882
        nd_set=1675632179.78
        
    ###L-band from Mel###
    if fname=='1630519596':
        nd_on_time=0.292376308037
        nd_cycle=19.4917538692
        nd_set=1630519596.16

    if fname=='1631379874':
        nd_on_time=0.292376308037
        nd_cycle=19.4917538692
        nd_set=1631379874.12
    
    if fname=='1631387336':
        nd_on_time=0.292376308037
        nd_cycle=19.4917538692
        nd_set=1631387336.24
    
    if fname=='1631659886':
        nd_on_time=0.584752616075
        nd_cycle=19.4917538692
        nd_set=1631659886.03
    
    if fname=='1631667564':
        nd_on_time=0.584752616075
        nd_cycle=19.4917538692
        nd_set=1631667564.78
    
    if fname=='1631552188':
        nd_on_time=0.584752616075
        nd_cycle=19.4917538692
        nd_set=1631552188.37
    
    if fname=='1631559762':
        nd_on_time=0.584752616075
        nd_cycle=19.4917538692
        nd_set=1631559762.04

    if fname=='1631724508':
        nd_on_time=0.584752616075
        nd_cycle=19.4917538692
        nd_set=1631724508.42
        
    if fname=='1631732038':
        nd_on_time=0.584752616075
        nd_cycle=19.4917538692
        nd_set=1631732038.16
    
    if fname=='1631810671':
        nd_on_time=0.584752616075
        nd_cycle=19.4917538692
        nd_set=1631810671.59
    
    if fname=='1631818149':
        nd_on_time=0.584752616075
        nd_cycle=19.4917538692
        nd_set=1631818149.54
    
    if fname=='1634835083':
        nd_on_time=0.584752616075
        nd_cycle=19.4917538692
        nd_set=1634835083.86
    
    if fname=='1631982988':
        nd_on_time=0.584752616075
        nd_cycle=19.4917538692
        nd_set=1631982988.79
        
    if fname=='1631990463':
        nd_on_time=0.584752616075
        nd_cycle=19.4917538692
        nd_set=1631990463.06
        
    if fname=='1632069690':
        nd_on_time=0.584752616075
        nd_cycle=19.4917538692
        nd_set=1632069690.79
        
    if fname=='1632077222':
        nd_on_time=0.584752616075
        nd_cycle=19.4917538692
        nd_set=1632077222.54
        
    if fname=='1632184922':
        nd_on_time=0.584752616075
        nd_cycle=19.4917538692
        nd_set=1632184922.75
        
    if fname=='1632505883':
        nd_on_time=0.584752616075
        nd_cycle=19.4917538692
        nd_set=1632505883.88
        
    if fname=='1632760885':
        nd_on_time=0.584752616075
        nd_cycle=19.4917538692
        nd_set=1632760885.0
        
    if fname=='1633365980':
        nd_on_time=0.584752616075
        nd_cycle=19.4917538692
        nd_set=1633365980.96
    
    if fname=='1633970780':
        nd_on_time=0.584752616075
        nd_cycle=19.4917538692
        nd_set=1633970780.76
        
    if fname=='1634252028':
        nd_on_time=0.584752616075
        nd_cycle=19.4917538692
        nd_set=1634252028.72
        
    if fname=='1634402485':
        nd_on_time=0.584752616075
        nd_cycle=19.4917538692
        nd_set=1634402485.21
        
    if fname=='1634748682':
        nd_on_time=0.584752616075
        nd_cycle=19.4917538692
        nd_set=1634748682.76
    
    if fname=='1637346562':
        nd_on_time=0.584752616075
        nd_cycle=19.4917538692
        nd_set=1637346562.68
        
    if fname=='1637354605':
        nd_on_time=0.584752616075
        nd_cycle=19.4917538692
        nd_set=1637354605.11
        
    if fname=='1637691677':
        nd_on_time=0.584752616075
        nd_cycle=19.4917538692
        nd_set=1637691677.34
        
    if fname=='1637699408':
        nd_on_time=0.584752616075
        nd_cycle=19.4917538692
        nd_set=1637699408.3
        
    if fname=='1638130295':
        nd_on_time=0.584752616075
        nd_cycle=19.4917538692
        nd_set=1638130295.17
        
    if fname=='1638294319':
        nd_on_time=0.584752616075
        nd_cycle=19.4917538692
        nd_set=1638294319.33
        
    if fname=='1638301944':
        nd_on_time=0.584752616075
        nd_cycle=19.4917538692
        nd_set=1638301944.04
        
    if fname=='1638386189':
        nd_on_time=0.584752616075
        nd_cycle=19.4917538692
        nd_set=1638386189.11
        
    if fname=='1638639082':
        nd_on_time=0.584752616075
        nd_cycle=19.4917538692
        nd_set=1638639082.19
        
    if fname=='1638647186':
        nd_on_time=0.584752616075
        nd_cycle=19.4917538692
        nd_set=1638647186.71
        
    if fname=='1638898468':
        nd_on_time=0.584752616075
        nd_cycle=19.4917538692
        nd_set=1638898468.77
        
    if fname=='1639157507':
        nd_on_time=0.584752616075
        nd_cycle=19.4917538692
        nd_set=1639157507.43
        
    if fname=='1639331184':
        nd_on_time=0.584752616075
        nd_cycle=19.4917538692
        nd_set=1639331184.4
        
    if fname=='1639935088':
        nd_on_time=0.584752616075
        nd_cycle=19.4917538692
        nd_set=1639935088.65
        
    if fname=='1640540184':
        nd_on_time=0.584752616075
        nd_cycle=19.4917538692
        nd_set=1640540184.68
        
    if fname=='1640712986':
        nd_on_time=0.584752616075
        nd_cycle=19.4917538692
        nd_set=1640712986.13
        
    if fname=='1640799689':
        nd_on_time=0.584752616075
        nd_cycle=19.4917538692
        nd_set=1640799689.6
    ###from Mel###
    
    if nd_on_time==-999:
        print ('# No record, can ask < astro.jywang@gmail.com >')
        
    return nd_on_time,nd_cycle,nd_set

def cal_nd_edges(timestamps,nd_set,nd_cycle,nd_on_time):
    nd_on_edge,nd_off_edge=[],[]
    for i in range(1000):
        edge1=nd_set+nd_cycle*i
        edge2=edge1+nd_on_time

        if edge1>timestamps[-1]:
            print ('edge number 0-'+str(i-1))
            break
        nd_on_edge.append(edge1)
        nd_off_edge.append(edge2)
    nd_on_edge,nd_off_edge=np.array(nd_on_edge),np.array(nd_off_edge)
    return nd_on_edge,nd_off_edge

def cal_nd_ratio(timestamps, nd_on_time, nd_on_edge, dump_period):
    nd_ratio=np.zeros_like(timestamps)
    nd_1x=[]

    for i in range(len(nd_on_edge)):
        on_edge_local=nd_on_edge[i]
        gap_list=abs(on_edge_local-timestamps)
        dp_gap_min=np.where(gap_list==np.min(gap_list))[0][0]
        gap_min=gap_list[dp_gap_min]
        #print i, dp_gap_min, gap_min,
        if gap_min > dump_period/2. :
            print ('*** diode '+str(i)+' was fired out of timestamps list: '+str(on_edge_local-timestamps[0])+' not in [0,'+str(timestamps[-1]-timestamps[0])+']')
        else:
            #time_edge=(timestamps[dp_gap_min]+timestamps[dp_gap_min+1])/2 
            time_edge=timestamps[dp_gap_min]+dump_period/2. 
            if time_edge - on_edge_local < nd_on_time: #diode in two dumps
                ratio1=(time_edge - on_edge_local)/dump_period
                nd_1x.append(dp_gap_min)
                nd_ratio[dp_gap_min]=ratio1
                if  dp_gap_min+1 <len(timestamps): #to make sure the second nd in the timestamps list    
                    ratio2=nd_on_time/dump_period-ratio1
                    nd_1x.append(dp_gap_min+1)
                    nd_ratio[dp_gap_min+1]=ratio2
                #print ratio1+ratio2,
                #print '=',
                #print ratio1,ratio2
            if time_edge - on_edge_local >= nd_on_time:  #diode in one dump
                ratio=nd_on_time/dump_period            
                nd_1x.append(dp_gap_min)            
                nd_ratio[dp_gap_min]=ratio
                #print ratio
                
    nd_0=np.where(nd_ratio==0)[0]
    #print nd_0
    assert(len(nd_0)+len(nd_1x)==len(timestamps))
    print ('# checked: len(nd_0)+len(nd_1x)==len(timestamps)')
    return nd_ratio, nd_0, nd_1x

def cal_nd_ratio_group(nd_ratio):
    nd_ratio_set=list(set(nd_ratio)) #check how many numbers for nd_ratio
    #print nd_ratio_set

    list_local=[]
    for i in nd_ratio_set:
        list_local.append(round(i,5)) #ignore the small differences
    nd_ratio_group=list(set(list_local))
    nd_ratio_group.sort()
    print (nd_ratio_group)
    return nd_ratio_group

