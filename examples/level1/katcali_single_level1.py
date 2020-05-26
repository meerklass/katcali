#!/usr/bin/env python2
import matplotlib
matplotlib.use('AGG')

#imports
import katdal
import numpy as np
import matplotlib.pylab as plt
import astropy.coordinates as ac
import functools
import healpy as hp
import optparse
import warnings
from matplotlib.backends.backend_pdf import PdfPages
import healpy as hp
from astropy import units as u
from matplotlib.offsetbox import AnchoredText
import time
import pickle
import sys
Tcmb=2.725
import katcali
import katcali.visualizer as kv
import katcali.models as km
import katcali.rfi as kr
import katcali.solver as ks
import katcali.io as kio
import katcali.label_dump as kl
import katcali.diode as kd
import katcali.filter as kf
import katcali.beam as kb
from astropy.coordinates import SkyCoord
from astropy import units as u

output_file='level1_output/'
# Select an observation block and load basic information in

fname=sys.argv[1]  

data=kio.load_data(fname)
print data

target,c0,bad_ants,flux_model=kio.check_ants(fname)
ch_plot=800
p_radec=np.loadtxt('radio_source.txt')

ant=sys.argv[2]

for pol in ['h','v']:
    print ant, pol
    ####################################################################
    #load data, labels, and parameters
    data.select(ants=ant,pol=pol)
    recv=ant+pol
    corr_id=kio.cal_corr_id(data,recv)
    assert(recv==data.corr_products[corr_id][0])
    assert(recv==data.corr_products[corr_id][1])
    print corr_id,recv
    vis,flags= kio.call_vis(fname,recv)
    vis_backup=vis.copy()
    ra,dec,az,el=kio.load_coordinates(data)
    timestamps,freqs=kio.load_tf(data)
    nd_set,nd_time,nd_cycle,nd_ratio=kio.load_ndparam(fname,data)
    dp_tt,dp_ss,dp_f,dp_t,dp_s=kl.label_dump_1ch(data,ant,pol,flags,ch_plot)
    dp_w=kl.select_waste(data,ant,pol)
    assert(np.shape(data)[2]==1)
    dp_sb=dp_ss[0]
    dp_se=dp_ss[-1]
    ang_deg=kio.load_ang_deg(ra,dec,c0)

    #read noise diode labels in 
    t_line=kd.cal_t_line(fname, timestamps,nd_set, nd_cycle, data.dump_period)
    #mark,nd_1_det,nd_1a_det,nd_1b_det,lmin,lmax=kd.label_nd_injection(fname,vis, timestamps, dp_ss, data.dump_period)
    nd_1a,nd_1b,nd_1,nd_0=kd.call_nd_1_list(fname,timestamps)
    nd_s1a,nd_s1b,nd_s1,nd_s0=kd.cal_nds_list(dp_ss,nd_1a,nd_1b,nd_1,nd_0)#dp_ss here, not dp_s
    nd_t1a,nd_t1b,nd_t1,nd_t0=kd.cal_ndt_list(dp_tt,nd_1a,nd_1b,nd_1,nd_0)#dp_tt here, not dp_t
    #nd_labels=nd_s1a,nd_s1b,nd_s1,nd_s0,nd_t1a,nd_t1b,nd_t1,nd_t0
    nd_label=dp_s,dp_t,nd_1a,nd_1b,nd_1,nd_0 ###only for rfi flagging

    #noise diode injection along time ***plot is very long***
    ch_plot0=800
    plt.figure(figsize=(9,4))
    plt.step(timestamps-timestamps[0],(vis[:,ch_plot0]),'gray',where='mid',lw=2)
    plt.plot(timestamps-timestamps[0],(vis[:,ch_plot0]),'b.')
    plt.plot(timestamps[nd_1a]-timestamps[0],(vis[nd_1a,ch_plot0]),'r.')
    plt.plot(timestamps[nd_1b]-timestamps[0],(vis[nd_1b,ch_plot0]),'g.')
    plt.xlabel('time (s)',fontsize=14)
    plt.legend(['dump edge','vis in nd_0','vis in nd_1a','vis in nd_1b'],fontsize=10,ncol=4)
    plt.ylabel('raw visibility',fontsize=14)
    plt.title('raw vis at '+str(round(freqs[ch_plot0]/1e9,3))+' GHz, obs'+str(fname)+', '+str(recv),fontsize=14)
    plt.xlim(1500,1660)
    #plt.ylim(230,410)
    #plt.xticks(np.arange(0,timestamps[-1]-timestamps[0]+1, 50))
    plt.grid()
    plt.savefig(output_file+fname+'_'+recv+'_his_part.png',bbox_inches='tight')
    #plt.show()

    ch_plot_list=[600,800,1000,2200,2400,2600,2800,3000]
    plt.figure(figsize=(24,20))
    for i in range(len(ch_plot_list)):
        ch_plot1=ch_plot_list[i]
        p_data=vis[nd_s0,ch_plot1]
        plt.subplot(5,2,i+1)
        plt.scatter(ra[nd_s0],dec[nd_s0], c=p_data, vmin=p_data.min(),vmax=p_data.max(), s=8)
        if i>len(ch_plot_list)-3:
            plt.xlabel('R.A.')
        plt.ylabel('Dec')
        plt.title(fname+'_'+recv+'_ch'+str(ch_plot1))
        #plt.plot(p_radec[:,0],p_radec[:,1],'mo')
    plt.savefig(output_file+fname+'_'+recv+'_scatter_map0.png')
    #plt.show()

    plt.figure(figsize=(24,20))
    plt.subplots_adjust(wspace =0)
    for i in range(len(ch_plot_list)):
        ch_plot1=ch_plot_list[i]
        p_data=vis[nd_s0,ch_plot1]
        plt.subplot(5,2,i+1)
        kv.plot_data(ra[nd_s0],dec[nd_s0], p_data,gsize=90)
        plt.plot(p_radec[:,0],p_radec[:,1],'mo')
        if i>len(ch_plot_list)-3:
            plt.xlabel('R.A.')
        plt.ylabel('Dec')
        plt.title(fname+'_'+recv+'_ch'+str(ch_plot1))
        #plt.plot(p_radec[:,0],p_radec[:,1],'mo')
    plt.savefig(output_file+fname+'_'+recv+'_map0.png')
    #plt.show()

    p = SkyCoord(data.ra*u.deg,  data.dec*u.deg, frame='icrs')
    dp_ptr_list=[]

    for i in range(len(p_radec)):
        #print i
        p_ra,p_dec=p_radec[i]
        c = SkyCoord(p_ra*u.deg,  p_dec*u.deg, frame='icrs')
        #print c 
        p_ang=(c.separation(p)/u.deg)[:,0]
        #print p_ang
        dp_l=np.where(p_ang<.5)[0]
        #print dp_l
        for j in range(len(dp_l)):
            if dp_l[j]>dp_sb and dp_l[j]<=dp_se:
                dp_ptr_list.append(dp_l[j])

    list(set(dp_ptr_list))
    dp_ptr_list.sort()
    dp_s0_ptr=list(set(dp_ptr_list).intersection(nd_s0))

    nd_t0_ca,nd_t0_cb,nd_t1a_ca,nd_t1a_cb,nd_t1b_ca,nd_t1b_cb=kd.cal_nd_t_c_list(nd_t0,nd_t1a, nd_t1b, dp_sb,dp_se)

    First_Thresholds=[] 
    for l in [nd_s0,nd_s1a,nd_s1b]:
        #print np.ma.mean(vis[l,:]), np.ma.median(vis[l,:])
        First_Thresholds.append(np.ma.median(vis[l,:])/10.)
    for l in [nd_t0_ca,nd_t1a_ca,nd_t1b_ca,nd_t0_cb,nd_t1a_cb,nd_t1b_cb]:
        #print np.ma.mean(vis[l,:]), np.ma.median(vis[l,:])
        First_Thresholds.append(np.ma.median(vis[l,:])/5.)

    print First_Thresholds

    vis_clean=kr.vis_flag(vis,flags,nd_label, dp_w, First_Thresholds,flag_step=1)

    plt.figure(figsize=(24,20))
    for i in range(len(ch_plot_list)):
        ch_plot1=ch_plot_list[i]
        p_data=vis_clean[nd_s0,ch_plot1]
        plt.subplot(5,2,i+1)
        plt.scatter(ra[nd_s0],dec[nd_s0], c=p_data, vmin=p_data.min(),vmax=p_data.max(), s=8)
        if i>len(ch_plot_list)-3:
            plt.xlabel('R.A.')
        plt.ylabel('Dec')
        plt.title(fname+'_'+recv+'_ch'+str(ch_plot1))
        #plt.plot(p_radec[:,0],p_radec[:,1],'mo')
    plt.savefig(output_file+fname+'_'+recv+'_scatter_map1.png')
    #plt.show()

    plt.figure(figsize=(24,20))
    plt.subplots_adjust(wspace =0)
    for i in range(len(ch_plot_list)):
        ch_plot1=ch_plot_list[i]
        p_data=vis_clean[nd_s0,ch_plot1]
        plt.subplot(5,2,i+1)
        kv.plot_data(ra[nd_s0],dec[nd_s0], p_data,gsize=90)
        plt.plot(p_radec[:,0],p_radec[:,1],'mo')
        if i>len(ch_plot_list)-3:
            plt.xlabel('R.A.')
        plt.ylabel('Dec')
        plt.title(fname+'_'+recv+'_ch'+str(ch_plot1))
        #plt.plot(p_radec[:,0],p_radec[:,1],'mo')
    plt.savefig(output_file+fname+'_'+recv+'_map1.png')
    #plt.show()

    ###second rfi
    from scipy.interpolate import Rbf
    rbf = Rbf(range(4096) ,np.ma.mean(vis[nd_0,:],axis=0), smooth=1)
    cur=rbf(range(4096))

    vis_clean2=vis_clean.copy()
    for ch_i in range(4096):
        vis_clean2[:,ch_i]=vis_clean2[:,ch_i]/cur[ch_i] #rescaled!!!
    print vis_clean2.mean()

    vis_clean2_part1=vis_clean2[:,500:1100]
    vis_clean2_part2=vis_clean2[:,2100:3150]

    vis_clean2_part1.mask[dp_ptr_list,:]=True
    vis_clean2_part2.mask[dp_ptr_list,:]=True

    Second_Thresholds1=[]
    for l in [nd_s0,nd_s1a,nd_s1b]:
        #print np.ma.mean(vis_clean2_part1[l,:]), np.ma.median(vis_clean2_part1[l,:])
        Second_Thresholds1.append(np.ma.median(vis_clean2_part1[l,:])/5.)
    for l in [nd_t0_ca,nd_t1a_ca,nd_t1b_ca,nd_t0_cb,nd_t1a_cb,nd_t1b_cb]:
        #print np.ma.mean(vis_clean2_part1[l,:]), np.ma.median(vis_clean2_part1[l,:])
        Second_Thresholds1.append(np.ma.median(vis_clean2_part1[l,:])/5.)

    Second_Thresholds2=[]
    for l in [nd_s0,nd_s1a,nd_s1b]:
        #print np.ma.mean(vis_clean2_part2[l,:]), np.ma.median(vis_clean2_part2[l,:])
        Second_Thresholds2.append(np.ma.median(vis_clean2_part2[l,:])/5.)
    for l in [nd_t0_ca,nd_t1a_ca,nd_t1b_ca,nd_t0_cb,nd_t1a_cb,nd_t1b_cb]:
        #print np.ma.mean(vis_clean2_part2[l,:]), np.ma.median(vis_clean2_part2[l,:])
        Second_Thresholds2.append(np.ma.median(vis_clean2_part2[l,:])/5.)

    print Second_Thresholds1
    print Second_Thresholds2

    #rfi flagging for raw vis data
    vis_clean2_part1=kr.vis_flag(vis_clean2_part1,vis_clean2_part1.mask,nd_label, dp_w, Second_Thresholds1,flag_step=2)
    vis_clean2_part2=kr.vis_flag(vis_clean2_part2,vis_clean2_part2.mask,nd_label, dp_w, Second_Thresholds2,flag_step=2)

    vis_clean2=vis_clean.copy() #reset vis_clean2
    vis_clean2.mask[:,500:1100]=vis_clean2_part1.mask
    vis_clean2.mask[:,2100:3150]=vis_clean2_part2.mask
    vis_clean2.mask[dp_ptr_list,:]=vis_clean.mask[dp_ptr_list,:]

    vis_clean2=kr.clean_bad_ratio(vis_clean2)

    plt.figure(figsize=(24,20))
    for i in range(len(ch_plot_list)):
        ch_plot1=ch_plot_list[i]
        p_data=vis_clean2[nd_s0,ch_plot1]
        plt.subplot(5,2,i+1)
        plt.scatter(ra[nd_s0],dec[nd_s0], c=p_data, vmin=p_data.min(),vmax=p_data.max(), s=8)
        if i>len(ch_plot_list)-3:
            plt.xlabel('R.A.')
        plt.ylabel('Dec')
        plt.title(fname+'_'+recv+'_ch'+str(ch_plot1))
        #plt.plot(p_radec[:,0],p_radec[:,1],'mo')
    plt.savefig(output_file+fname+'_'+recv+'_scatter_map2.png')
    #plt.show()

    plt.figure(figsize=(24,20))
    plt.subplots_adjust(wspace =0)
    for i in range(len(ch_plot_list)):
        ch_plot1=ch_plot_list[i]
        p_data=vis_clean2[nd_s0,ch_plot1]
        plt.subplot(5,2,i+1)
        kv.plot_data(ra[nd_s0],dec[nd_s0], p_data,gsize=90)
        plt.plot(p_radec[:,0],p_radec[:,1],'mo')
        if i>len(ch_plot_list)-3:
            plt.xlabel('R.A.')
        plt.ylabel('Dec')
        plt.title(fname+'_'+recv+'_ch'+str(ch_plot1))
        #plt.plot(p_radec[:,0],p_radec[:,1],'mo')
    plt.savefig(output_file+fname+'_'+recv+'_map2.png')
    #plt.show()

    vis_clean=vis_clean2.copy() #update vis_clean

    #compare vis before and after rfi flagging
    plt.figure(figsize=(14,5.4))
    plt.subplot(121)
    plt.imshow(vis_backup,aspect='auto')
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.ylabel('time',fontsize=12)
    plt.xlabel('channel',fontsize=12)
    plt.title('raw vis of '+str(fname)+', '+str(recv),fontsize=16, y=1.12)
    plt.colorbar()
    plt.twiny()
    plt.xticks(fontsize=12)
    plt.imshow(vis_backup,aspect='auto',extent=(data.freqs[0]/1e9,data.freqs[-1]/1e9,len(timestamps),0))
    plt.xlabel('Freq (GHz)',fontsize=12)
    plt.subplot(122)
    plt.imshow(vis_clean,aspect='auto')
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.ylabel('time',fontsize=12)
    plt.xlabel('channel',fontsize=12)
    plt.title('clean vis of '+str(fname)+', '+str(recv),fontsize=16, y=1.12)
    plt.colorbar()
    plt.twiny()
    plt.xticks(fontsize=12)
    plt.imshow(vis_clean,aspect='auto',extent=(data.freqs[0]/1e9,data.freqs[-1]/1e9,len(timestamps),0))
    plt.xlabel('Freq (GHz)',fontsize=12)
    #plt.savefig(str(fname)+'raw_vis.pdf', bbox_inches='tight')
    plt.savefig(output_file+str(fname)+'_'+recv+'_raw_vis.png', bbox_inches='tight')
    #plt.show()

    plt.figure(figsize=(12,5))
    plt.subplot(121)
    plt.imshow(vis_clean,aspect='auto')
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.ylabel('time',fontsize=12)
    plt.xlabel('channel',fontsize=12)
    plt.xlim(550,1050)
    plt.colorbar()
    plt.subplot(122)
    plt.imshow(vis_clean,aspect='auto')
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.ylabel('time',fontsize=12)
    plt.xlabel('channel',fontsize=12)
    plt.xlim(2150,3100)
    plt.colorbar()
    plt.savefig(output_file+str(fname)+'_'+recv+'_raw_vis_cut.png', bbox_inches='tight')
    #plt.show()

    d={}
    d['mask']=vis_clean.mask
    fs=open(output_file+str(fname)+'_'+str(recv)+'_mask','wb')
    pickle.dump(d,fs,protocol=2)
    fs.close()

    if pol=='h':
        vis_clean_h=vis_clean.copy()
    if pol=='v':
        vis_clean_v=vis_clean.copy()

print ant
vis_div=vis_clean_h/vis_clean_v #intersection mask
vis_div=kr.clean_bad_ratio(vis_div,ratio_t=0.6,ratio_ch=0.5) #clean bad ratio part

vis_clean_hh=np.ma.array(vis_clean_h,mask=vis_div.mask)
vis_clean_vv=np.ma.array(vis_clean_v,mask=vis_div.mask)

d={}
d['mask']=vis_div.mask
fs=open(output_file+str(fname)+'_'+str(ant)+'_mask','wb')
pickle.dump(d,fs,protocol=2)
fs.close()


plt.figure(figsize=(14,5.4))
plt.subplot(121)
plt.imshow(vis_clean_hh,aspect='auto')
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.ylabel('time',fontsize=12)
plt.xlabel('channel',fontsize=12)
plt.title('clean vis of '+str(fname)+', '+str(ant)+'h',fontsize=16, y=1.12)
plt.colorbar()
plt.twiny()
plt.xticks(fontsize=12)
plt.imshow(vis_clean_hh,aspect='auto',extent=(data.freqs[0]/1e9,data.freqs[-1]/1e9,len(timestamps),0))
plt.xlabel('Freq (GHz)',fontsize=12)
plt.subplot(122)
plt.imshow(vis_clean_vv,aspect='auto')
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.ylabel('time',fontsize=12)
plt.xlabel('channel',fontsize=12)
plt.title('clean vis of '+str(fname)+', '+str(ant)+'v',fontsize=16, y=1.12)
plt.colorbar()
plt.twiny()
plt.xticks(fontsize=12)
plt.imshow(vis_clean_vv,aspect='auto',extent=(data.freqs[0]/1e9,data.freqs[-1]/1e9,len(timestamps),0))
plt.xlabel('Freq (GHz)',fontsize=12)
#plt.savefig(str(fname)+'raw_vis.pdf', bbox_inches='tight')
plt.savefig(output_file+str(fname)+'_'+ant+'_clean_vis.png', bbox_inches='tight')
#plt.show()

for pol in ['h','v']:
    recv=ant+pol
    if pol=='h':
        vis_clean_pol=vis_clean_hh
    if pol=='v':
        vis_clean_pol=vis_clean_vv
        
    plt.figure(figsize=(24,20))
    for i in range(len(ch_plot_list)):
        ch_plot1=ch_plot_list[i]
        p_data=vis_clean_pol[nd_s0,ch_plot1]
        plt.subplot(5,2,i+1)
        plt.scatter(ra[nd_s0],dec[nd_s0], c=p_data, vmin=p_data.min(),vmax=p_data.max(), s=8)
        if i>len(ch_plot_list)-3:
            plt.xlabel('R.A.')
        plt.ylabel('Dec')
        plt.title(fname+'_'+recv+'_ch'+str(ch_plot1))
        #plt.plot(p_radec[:,0],p_radec[:,1],'mo')
    plt.savefig(output_file+fname+'_'+recv+'_scatter_map.png')
    #plt.show()

    plt.figure(figsize=(24,20))
    plt.subplots_adjust(wspace =0)
    for i in range(len(ch_plot_list)):
        ch_plot1=ch_plot_list[i]
        p_data=vis_clean_pol[nd_s0,ch_plot1]
        plt.subplot(5,2,i+1)
        kv.plot_data(ra[nd_s0],dec[nd_s0], p_data,gsize=90)
        plt.plot(p_radec[:,0],p_radec[:,1],'mo')
        if i>len(ch_plot_list)-3:
            plt.xlabel('R.A.')
        plt.ylabel('Dec')
        plt.title(fname+'_'+recv+'_ch'+str(ch_plot1))
        #plt.plot(p_radec[:,0],p_radec[:,1],'mo')
    plt.savefig(output_file+fname+'_'+recv+'_map.png')
    #plt.show()

################################################
#### Jingying Wang <astro.jywang@gmail.com> ####
################################################
