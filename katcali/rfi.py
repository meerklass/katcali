from . import SEEK_sum_threshold as sum_threshold
import numpy as np
import matplotlib.pylab as plt
import numpy.ma as ma
from . import filter as kf
from . import diode as kd 
from . import label_dump as kl

#param set ref: https://seek.readthedocs.io/en/latest/_modules/seek/mitigation/sum_threshold.html
def seek_rfi_mask(autodata, First_Threshold, sm_kwargs_para=(80,80,40,40),di_kwargs_para=(25,30)):
    rfi_mask = sum_threshold.get_rfi_mask(tod=autodata.astype('float'),
                                              mask=autodata.mask.astype('bool'),
                                              chi_1=First_Threshold,
                                              #sm_kwargs=sum_threshold.get_sm_kwargs(40, 20, 15, 7.5),
                                              sm_kwargs=sum_threshold.get_sm_kwargs(*sm_kwargs_para),
                                              #di_kwargs=sum_threshold.get_di_kwrags(3, 7),
                                              di_kwargs=sum_threshold.get_di_kwrags(*di_kwargs_para),
                                              plotting=False)

    tod_mask = ma.masked_array(autodata, mask=rfi_mask) 
    filval = -99999
    filval_off = np.std(tod_mask.compressed())
    filval_off2 = np.std(tod_mask)/np.mean(tod_mask)
    print ('------------------------------------------')
    print ('sm_kwargs_para: '+str(sm_kwargs_para))
    print ('di_kwargs_para: '+str(di_kwargs_para))
    #print 'Std of the non-masked elements of tod is ' + str(round(filval_off,2))
    print ('Std/Mean of the non-masked elements of tod is ' + str(round(filval_off2,2)))
    
    return tod_mask

###########################################full vis mask###############################################################
def vis_flag(vis_backup,flags,nd_label0, dp_w, First_Thresholds,flag_step, key=0):
    dp_s,dp_t,nd_1a,nd_1b,nd_1,nd_0=nd_label0
    nd_s1a,nd_s1b,nd_s1,nd_s0=kd.cal_nds_list(dp_s,nd_1a,nd_1b,nd_1,nd_0)#dp_s here, not dp_ss
    nd_t1a,nd_t1b,nd_t1,nd_t0=kd.cal_ndt_list(dp_t,nd_1a,nd_1b,nd_1,nd_0)#dp_t here, not dp_tt

    
    nd_t0_ca,nd_t0_cb,nd_t1a_ca,nd_t1a_cb,nd_t1b_ca,nd_t1b_cb=kd.cal_nd_t_c_list(nd_t0,nd_t1a, nd_t1b, dp_s[0],dp_s[-1])
        
    print ('### load flags ###')
    vis=np.ma.array(vis_backup,mask=flags)
    n=np.shape(np.where(flags==True))[1]
    print ('data.flags ratio:')
    print (float(n)/np.shape(vis)[0]/np.shape(vis)[1])

    print ('###mask data not track/scan###') 
    for i in range(np.shape(vis)[0]):
        if i in dp_w:
            vis[i,:].mask='True'

         
    print ('###SEEK flagging###')
    print ('First_Thresholds for nd_s0, nd_s1a, nd_s1b, nd_t0_ca, nd_t1a_ca, nd_t1b_ca, nd_t0_cb, nd_t1a_cb,nd_t1b_cb= '+ str(First_Thresholds))
    print ('flag_step: '+str(flag_step))
    
    if flag_step==1:
        sm_kwargs_para_nd0=(60, 50, 30, 25)
        di_kwargs_para_nd0=(20, 20)
        sm_kwargs_para_nd1=(60, 50, 30, 25)
        di_kwargs_para_nd1=(1, 20)
    if flag_step==2:
        sm_kwargs_para_nd0=(40, 20, 20, 10)
        di_kwargs_para_nd0=(5, 10)
        sm_kwargs_para_nd1=(40, 20, 20, 10)
        di_kwargs_para_nd1=(1, 10)
    
    if len(nd_s0)>0:
        print ('#flagging scan and nd_0')
        vis_s0_m=seek_rfi_mask(vis[nd_s0,:], First_Thresholds[0], sm_kwargs_para=sm_kwargs_para_nd0,di_kwargs_para=di_kwargs_para_nd0)
    if len(nd_s1a)>0:
        print ('#flagging scan and nd_1a')
        vis_s1a_m=seek_rfi_mask(vis[nd_s1a,:], First_Thresholds[1], sm_kwargs_para=sm_kwargs_para_nd1,di_kwargs_para=di_kwargs_para_nd1)
    if len(nd_s1b)>0:    
        print ('#flagging scan and nd_1b')
        vis_s1b_m=seek_rfi_mask(vis[nd_s1b,:], First_Thresholds[2], sm_kwargs_para=sm_kwargs_para_nd1,di_kwargs_para=di_kwargs_para_nd1)
    if len(nd_t0_ca)>0:
        print ('#flagging track, nd_0, and before scan')
        vis_t0_ca_m=seek_rfi_mask(vis[nd_t0_ca,:], First_Thresholds[3], sm_kwargs_para=sm_kwargs_para_nd0,di_kwargs_para=di_kwargs_para_nd0)
    if len(nd_t1a_ca)>0:    
        print ('#flagging track, nd_1a, and before scan')
        vis_t1a_ca_m=seek_rfi_mask(vis[nd_t1a_ca,:], First_Thresholds[4], sm_kwargs_para=sm_kwargs_para_nd1,di_kwargs_para=di_kwargs_para_nd1)
    if len(nd_t1b_ca)>0:
        print ('#flagging track, nd_1b, and before scan')
        vis_t1b_ca_m=seek_rfi_mask(vis[nd_t1b_ca,:], First_Thresholds[5], sm_kwargs_para=sm_kwargs_para_nd1,di_kwargs_para=di_kwargs_para_nd1)
    if len(nd_t0_cb)>0:
        print ('#flagging track, nd_0, and after scan')
        vis_t0_cb_m=seek_rfi_mask(vis[nd_t0_cb,:], First_Thresholds[6], sm_kwargs_para=sm_kwargs_para_nd0,di_kwargs_para=di_kwargs_para_nd0)
    if len(nd_t1a_cb)>0:
        print ('#flagging track, nd_1a, and after scan')
        vis_t1a_cb_m=seek_rfi_mask(vis[nd_t1a_cb,:], First_Thresholds[7], sm_kwargs_para=sm_kwargs_para_nd1,di_kwargs_para=di_kwargs_para_nd1)
    if len(nd_t1b_cb)>0:
        print ('#flagging track, nd_1b, and after scan')
        vis_t1b_cb_m=seek_rfi_mask(vis[nd_t1b_cb,:], First_Thresholds[8], sm_kwargs_para=sm_kwargs_para_nd1,di_kwargs_para=di_kwargs_para_nd1)
    print ('---------------------------------------------------')        
    print ('#put all parts together')
    vis_clean=np.ma.array(np.zeros_like(vis),mask=True)
    ####scan part
    vis_clean[nd_s0,:]=vis_s0_m
    if len(nd_s1a)>0:
        vis_clean[nd_s1a,:]=vis_s1a_m
    if len(nd_s1b)>0:
        vis_clean[nd_s1b,:]=vis_s1b_m
    ####track part a
    vis_clean[nd_t0_ca,:]=vis_t0_ca_m
    if len(nd_t1a_ca)>0:
        vis_clean[nd_t1a_ca,:]=vis_t1a_ca_m
    if len(nd_t1b_ca)>0:
        vis_clean[nd_t1b_ca,:]=vis_t1b_ca_m
    ####track part b
    vis_clean[nd_t0_cb,:]=vis_t0_cb_m
    if len(nd_t1a_cb)>0:
        vis_clean[nd_t1a_cb,:]=vis_t1a_cb_m
    if len(nd_t1b_cb)>0:
        vis_clean[nd_t1b_cb,:]=vis_t1b_cb_m
    
    print ('#checking neighbours')  
    #if the neighbor diode off is masked, the diode on should be masked also
    for ch_i in range(np.shape(vis)[1]):
        for i in nd_s1a:
            if vis_clean.mask[i-1,ch_i]==True:
                vis_clean.mask[i,ch_i]=True

        for i in nd_s1b:
            if vis_clean.mask[i+1,ch_i]==True:
                vis_clean.mask[i,ch_i]=True

        for i in nd_t1a:
            if vis_clean.mask[i-1,ch_i]==True:
                vis_clean.mask[i,ch_i]=True

        for i in nd_t1b:
            if vis_clean.mask[i+1,ch_i]==True:
                vis_clean.mask[i,ch_i]=True


    print ('#cleaning the bad ratio part')
    if key==0:
        vis_clean2=clean_bad_ratio(vis_clean)
    if key==-1:
        print ('## ratio clean is cancelled!!!')
        vis_clean2=clean_bad_ratio(vis_clean, ratio_t=1.,ratio_ch=1.)
    return vis_clean2
###################################################################################
def clean_bad_ratio(vis,ratio_t=0.4,ratio_ch=0.5):
    vis2=vis.copy()
    t_len=np.shape(vis)[0]
    ch_len=np.shape(vis)[1]

    for i in range(ch_len):
        ratio=float(np.array(vis.mask[:,i]==True).sum())/t_len
        if ratio>ratio_ch:
            vis2.mask[:,i]=True
    for i in range(t_len):
        ratio=float(np.array(vis.mask[i,:]==True).sum())/ch_len
        if ratio>ratio_t: 
            vis2.mask[i,:]=True
    return vis2

##############single channel mask#######################################
def vis_flag_1ch(vis_clean,nd_s0,nd_s1a,nd_s1b,ch, sigma=5):
    
    print ('group shape (with flags):')
    print (len(nd_s0), len(nd_s1a),len(nd_s1b))

    nd_s0_clean=[]
    for i in nd_s0:
        if vis_clean.mask[i,ch]==False and vis_clean[i,ch]>0:
            nd_s0_clean.append(i)

    for iter in range(5):
        print ('iter='+str(iter))
    
        if iter==0:
            nd_s0_i=nd_s0_clean
            l1=len(nd_s0_clean)
        if iter>0:
            nd_s0_i=nd_s0_clean2
            l1=len(nd_s0_clean2)
            
        nd_s0_clean2=kf.curve_filter(nd_s0_i,vis_clean[nd_s0_i,ch],sigma=sigma) ###output is the data want to keep
        l2=len(nd_s0_clean2) #end len
        #print l1,l2
        
        if l1==l2:
            break
    nd_s0_clean=nd_s0_clean2
    
    #update vis_clean
    for i in nd_s0:
        if i not in nd_s0_clean:
            #print i
            vis_clean.mask[i,ch]=True 
            
    for i in nd_s1a:
        if vis_clean.mask[i-1,ch]==True:
            vis_clean.mask[i,ch]=True
                
    for i in nd_s1b:
        if vis_clean.mask[i+1,ch]==True:
            vis_clean.mask[i,ch]=True
    return vis_clean,nd_s0_clean


##########################added for 2021 data#######################################
def vis_flag_v2(data,flags,ch_ref,ant,pol,vis,timestamps, nd_on_time, nd_on_edge, dump_period, ang_deg, flag_step, Threshold_factor1, Threshold_factor2, ratio_clean_key=0, plt_key=0):
    
    dp_tt,dp_ss,dp_f,dp_w, dp_t,dp_s,dp_slew,dp_stop=kl.cal_dp_label(data,flags,ant,pol,ch_ref,ang_deg)
    nd_ratio,nd_0, nd_1x=kd.cal_nd_ratio(timestamps, nd_on_time, nd_on_edge, dump_period)
    nd_t0,nd_t1x,nd_s0,nd_s1x,nd_t0_ca,nd_t0_cb,nd_t1x_ca,nd_t1x_cb=kl.cal_label_intersec(dp_tt,dp_ss,nd_0,nd_1x)
    labels_1x=kl.cal_label_intersec_complex(dp_tt,dp_ss,nd_0,nd_1x,nd_ratio)
    
    labels=labels_1x+[list(nd_t0_ca)]+[list(nd_s0)]+[list(nd_t0_cb)]
    group_num=len(labels)
    group_num_1x=len(labels_1x)
    print ('group number in total: '+str(group_num))
    
    print ('### load flags ###')
    if np.shape(flags)==np.shape(vis):
        print ('### RFI flagging for full vis ### ')
        vis=np.ma.array(vis,mask=flags)
        n=np.shape(np.where(flags==True))[1]
        print ('data.flags ratio:')
        print (float(n)/np.shape(vis)[0]/np.shape(vis)[1])
    if np.shape(flags)> np.shape(vis):
        print ('### RFI flagging for small region ### ')

    print ('### mask data not track/scan ###') 
    for i in range(np.shape(vis)[0]):
        if i in dp_w:
            vis[i,:].mask='True'

         
    print ('### SEEK flagging ###')
    print ('flag_step: '+str(flag_step))
   
    #nd_ratio_group=cal_nd_ratio_group(nd_ratio)


    vis_clean=np.ma.array(np.zeros_like(vis),mask=True)
    dump_checker=[]
    
    for l_i in range(group_num):
        
        l=labels[l_i]
        print ('group '+str(l_i)+', dp num '+str(len(l)))
        if len(l)>0 and (vis[l,:].mask==True).all()==False:
            
            print ('data shape is '+str(np.shape(vis[l,:])))
            
            if len(l)>3:
                if flag_step==1:
                    if l_i<group_num_1x:
                        sm_kwargs_para_nd=(60, 50, 30, 25)
                        di_kwargs_para_nd=(1, 20)
                        Threshold_factor=Threshold_factor2
                        print ('Threshold_factor2='+str(Threshold_factor2)+' applied on nd_1x')
                    if l_i==group_num_1x+1:               
                        sm_kwargs_para_nd=(60, 50, 30, 25)
                        di_kwargs_para_nd=(20, 20)
                        Threshold_factor=Threshold_factor1
                        print ('Threshold_factor1='+str(Threshold_factor1)+' applied on nd_s0')
                    if l_i==group_num_1x or l_i==group_num_1x+2:               
                        sm_kwargs_para_nd=(60, 50, 30, 25)
                        di_kwargs_para_nd=(20, 20)
                        Threshold_factor=Threshold_factor1/2.
                        print ('Threshold_factor1='+str(Threshold_factor1/2.)+' applied on nd_t0')

                if flag_step==2:
                    if l_i<group_num_1x:
                        sm_kwargs_para_nd=(40, 20, 20, 10)
                        di_kwargs_para_nd=(1, 10)
                        Threshold_factor=Threshold_factor2
                        print ('Threshold_factor2='+str(Threshold_factor2)+' applied on nd_1x')
                    if l_i==group_num_1x+1:  
                        sm_kwargs_para_nd=(40, 20, 20, 10)
                        di_kwargs_para_nd=(5, 10)
                        Threshold_factor=Threshold_factor1
                        print ('Threshold_factor1='+str(Threshold_factor1)+' applied on nd_s0')
                    if l_i==group_num_1x or l_i==group_num_1x+2:  
                        sm_kwargs_para_nd=(40, 20, 20, 10)
                        di_kwargs_para_nd=(5, 10)
                        Threshold_factor=Threshold_factor1/2.
                        print ('Threshold_factor1='+str(Threshold_factor1/2.)+' applied on nd_t0')


                Threshold_local=np.ma.median(vis[l,:])/Threshold_factor

                print ('Threshold_local is '+str(Threshold_local))


                print ('RFI flagging applied on group'+str(l_i))
                vis_m=seek_rfi_mask(vis[l,:], Threshold_local, sm_kwargs_para=sm_kwargs_para_nd, di_kwargs_para=di_kwargs_para_nd)
                vis_clean[l,:]=vis_m.copy()
                dump_checker.append(list(l))
                
                if plt_key==0:

                    plt.figure(figsize=(15,3))
                    plt.subplot(131)
                    plt.imshow(vis[l,:],aspect='auto')
                    plt.colorbar()
                    plt.title('before')
                    plt.subplot(132)
                    plt.imshow(vis_m,aspect='auto')
                    plt.colorbar()
                    plt.title('after')
                    plt.subplot(133)
                    plt.imshow(vis_clean,aspect='auto')
                    plt.colorbar()
                    plt.title('add into vis_clean')
                    plt.show()
                
            if len(l)<=3:
                vis_clean[l,:]=vis[l,:].copy()
                dump_checker.append(list(l))
                print ('*** few data, no RFI flagging applied ***')
                                           
        else:
            print ('*** no data, skipped ***')
            
    vis_clean2=vis_clean.copy()    
  
    print ('#checking neighbours')  
    #if the neighbor diode off is masked, the diode on should be masked also
    for ch_i in range(np.shape(vis_clean)[1]):
        for i in nd_1x:
            if i==0 and vis_clean.mask[i+1,ch_i]==True:
                vis_clean2.mask[i,ch_i]=True
            if i==np.shape(vis_clean)[0]-1 and vis_clean.mask[i-1,ch_i]==True:
                vis_clean2.mask[i,ch_i]=True
            if i>0 and i<np.shape(vis_clean)[0]-1: 
                if vis_clean.mask[i-1,ch_i]==True or vis_clean.mask[i+1,ch_i]==True:
                    vis_clean2.mask[i,ch_i]=True
                    
    print ('#cleaning the bad ratio part')
    if ratio_clean_key==0:
        vis_clean3=clean_bad_ratio(vis_clean2)
    if ratio_clean_key==-1:
        print ('## ratio clean is cancelled!!!')
        vis_clean3=clean_bad_ratio(vis_clean2, ratio_t=1.,ratio_ch=1.)
    
    a=list(dump_checker).sort()
    b=(list(dp_tt)+list(dp_ss)).sort()
    assert(a==b)
    print ('# all dp_tt and dp_ss checked')
    
    return vis_clean3
