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
from astropy.coordinates import SkyCoord
from astropy import units as u
from scipy.interpolate import Rbf
print katcali.__version__

print 'start @ ' + time.asctime(time.localtime(time.time())) +'#'

print  plt.rcParams['font.size'], plt.rcParams[u'axes.linewidth'],plt.rcParams['lines.linewidth']
plt.rcParams['font.size'], plt.rcParams[u'axes.linewidth'],plt.rcParams['lines.linewidth'] = 14, 1.5, 1.5
#plt.rcParams['font.size'], plt.rcParams[u'axes.linewidth'],plt.rcParams['lines.linewidth'] = 10.0, 0.8, 1.5
print  plt.rcParams['font.size'], plt.rcParams[u'axes.linewidth'],plt.rcParams['lines.linewidth']

output_file="/scratch3/users/jywang/MeerKLASS2021/"

fname=sys.argv[1]  

data=kio.load_data(fname)
#print data
#data.obs_script_log

target,c0,bad_ants,flux_model=kio.check_ants(fname)

ants_good=[]
for i in np.array(kio.ant_list(data)):
    if i not in bad_ants:
        ants_good.append(i)
    else:
        print str(i) + ' is bad'
print fname
print ants_good

nd_on_time,nd_cycle,nd_set=kd.cal_nd_basic_para(fname)
print nd_on_time,nd_cycle,nd_set

ch_ref,ch_plot=800,800
p_radec=np.loadtxt('radio_source2021.txt')
ra_a4059,dec_a4059=-0.74042,-34.76056
ang_lim=.5

# Select ant and polarization, then load data in 

#select ant, polarization, and one channel to show data calibration
ant=sys.argv[2]
if ant in ants_good: 
    try:
        for pol in ['h','v']:
            data.select(ants=ant,pol=pol)
            recv=ant+pol
            corr_id=kio.cal_corr_id(data,recv)
            assert(recv==data.corr_products[corr_id][0])
            assert(recv==data.corr_products[corr_id][1])
            print corr_id,recv

            ra,dec,az,el=kio.load_coordinates(data)
            timestamps,freqs=kio.load_tf(data)
            dump_period=data.dump_period
            ang_deg=kio.load_ang_deg(ra,dec,c0)

            vis,flags= kio.call_vis(fname,recv)
            vis_backup=vis.copy()
            vis=np.ma.array(vis_backup,mask=flags)
            dp_tt,dp_ss,dp_f,dp_w, dp_t,dp_s,dp_slew,dp_stop=kl.cal_dp_label(data,flags,ant,pol,ch_ref,ang_deg)

            ra1=ra.copy()
            for i in range(len(ra)):
                if ra[i]>180:
                    ra1[i]=ra[i]-360

            print np.mean(ra),np.mean(ra1)
            ra=ra1
            print np.mean(ra),np.mean(ra1)

            nd_on_edge,nd_off_edge=kd.cal_nd_edges(timestamps,nd_set,nd_cycle,nd_on_time)
            print len(nd_on_edge),len(nd_off_edge)

            nd_ratio,nd_0, nd_1x=kd.cal_nd_ratio(timestamps, nd_on_time, nd_on_edge, dump_period)
            nd_t0,nd_t1x,nd_s0,nd_s1x,nd_t0_ca,nd_t0_cb,nd_t1x_ca,nd_t1x_cb=kl.cal_label_intersec(dp_tt,dp_ss,nd_0,nd_1x)

            labels_1x=kl.cal_label_intersec_complex(dp_tt,dp_ss,nd_0,nd_1x,nd_ratio)
            dp_sb=dp_ss[0]
            dp_se=dp_ss[-1]


            p = SkyCoord(data.ra*u.deg,  data.dec*u.deg, frame='icrs')
            dp_ptr_list=kl.cal_ptr_mask(p,p_radec,nd_s0, dp_sb,dp_se,ang_lim)

            # RFI flagging
            ## Basic RFI flagging (all channels)

            flag_step=1
            Threshold_factor1, Threshold_factor2=6,4 #diode off, on

            vis_clean=kr.vis_flag_v2(data, flags, ch_ref, ant,pol,vis_backup,timestamps, nd_on_time, nd_on_edge, dump_period, ang_deg, flag_step, Threshold_factor1, Threshold_factor2, ratio_clean_key=0, plt_key=-1)

            ###second rfi
            rbf = Rbf(range(4096) ,np.ma.mean(vis[nd_0,:],axis=0), smooth=1) #vis has no mask, flags hasn't applied.
            cur=rbf(range(4096))
            vis_clean2=vis_clean.copy()
            for ch_i in range(4096):
                vis_clean2[:,ch_i]=vis_clean2[:,ch_i]/cur[ch_i] #rescaled!!!
            print vis_clean2.mean()
            vis_clean2_part1=vis_clean2[:,500:1100]
            vis_clean2_part2=vis_clean2[:,2100:3150]

            vis_clean2_part1.mask[dp_ptr_list,:]=True
            vis_clean2_part2.mask[dp_ptr_list,:]=True

            flag_step=2
            Threshold_factor11,Threshold_factor22=1.1,1. #diode off, on
            #rfi flagging for raw vis data
            vis_clean2_part1=kr.vis_flag_v2(data, flags, ch_ref, ant,pol,vis_clean2_part1, timestamps, nd_on_time, nd_on_edge, dump_period, ang_deg, flag_step, Threshold_factor11, Threshold_factor22, ratio_clean_key=0, plt_key=-1)
            vis_clean2_part2=kr.vis_flag_v2(data, flags, ch_ref, ant,pol,vis_clean2_part2, timestamps, nd_on_time, nd_on_edge, dump_period, ang_deg, flag_step, Threshold_factor11, Threshold_factor22, ratio_clean_key=0, plt_key=-1)

            vis_clean2=vis_clean.copy() #reset vis_clean2
            vis_clean2.mask[:,500:1100]=vis_clean2_part1.mask
            vis_clean2.mask[:,2100:3150]=vis_clean2_part2.mask
            vis_clean2.mask[dp_ptr_list,:]=vis_clean.mask[dp_ptr_list,:]
            vis_clean2=kr.clean_bad_ratio(vis_clean2)
            vis_clean=vis_clean2.copy() #update vis_clean

            plot_gsize=60
            #raw map in different resolution
            plt.figure(figsize=(24,9))
            plt.subplots_adjust(wspace=0.,hspace=.2)
            plt.subplot(221)
            p_data=np.log10(vis_backup[nd_s0,ch_plot])
            plt.scatter(ra[nd_s0],dec[nd_s0], c=p_data, s=8)
            plt.title('Logged raw signal map of '+str(fname)+', '+str(recv)+ ' @'+ str(round(freqs[ch_plot]/1e6,0))+' MHz')
            #plt.xlabel('R.A. (J2000)')
            plt.ylabel('Dec (J2000)')
            clb = plt.colorbar()
            plt.subplot(222)
            p_data=vis_clean[nd_s0,ch_plot]
            plt.scatter(ra[nd_s0],dec[nd_s0], c=p_data, s=8)
            plt.title('RFI flagged raw signal map @'+ str(round(freqs[ch_plot]/1e6,0))+' MHz')
            #plt.xlabel('R.A. (J2000)')
            #plt.ylabel('Dec (J2000)')
            clb = plt.colorbar()
            plt.subplot(223)
            p_data=np.ma.array(np.log10(vis_backup[nd_s0,ch_plot]),mask=False) #for plot_mdata mask
            kv.plot_mdata(ra[nd_s0],dec[nd_s0], p_data,gsize=plot_gsize,  grid_method='linear', levels=6, x_mask=1, y_mask=2)
            #plt.plot(p_radec[:,0],p_radec[:,1],'mo')
            plt.title('smoothed map of above')
            plt.xlabel('R.A. (J2000)')
            plt.ylabel('Dec (J2000)')
            plt.subplot(224)
            p_data=vis_clean[nd_s0,ch_plot]
            kv.plot_mdata(ra[nd_s0],dec[nd_s0], p_data,gsize=plot_gsize,  grid_method='linear', levels=6, x_mask=1, y_mask=2)
            #plt.plot(p_radec[:,0],p_radec[:,1],'mo')
            plt.title('smoothed map of above')
            plt.xlabel('R.A. (J2000)')
            #plt.ylabel('Dec (J2000)')
            plt.savefig(output_file+'F_'+fname+'_'+recv+'_rfi_comp_map.png', bbox_inches='tight')
            #plt.show()

            d={}
            d['mask']=vis_clean.mask
            fs=open(output_file+fname+'_'+str(recv)+'_mask','wb')
            pickle.dump(d,fs,protocol=2)
            fs.close()


        #get vis
        recv1=ant+'h'
        recv2=ant+'v'
        vis_h,flags_h= kio.call_vis(fname,recv1)
        vis_v,flags_v= kio.call_vis(fname,recv2)

        #get mask
        d1 = pickle.load(open(output_file+fname+'_'+ant+'h_mask'))
        mask_h=d1['mask']
        d2 = pickle.load(open(output_file+fname+'_'+ant+'v_mask'))
        mask_v=d2['mask']
        #vis_clean
        vis_clean_h=np.ma.array(vis_h,mask=mask_h)
        vis_clean_v=np.ma.array(vis_v,mask=mask_v)
        print np.ma.mean(vis_clean_h),np.ma.mean(vis_clean_v)
        #get intersection mask
        vis_div=vis_clean_h/vis_clean_v #intersection mask
        vis_div=kr.clean_bad_ratio(vis_div,ratio_t=0.6,ratio_ch=0.5) #clean bad ratio part
        #update vis_clean
        vis_clean_hh=np.ma.array(vis_clean_h,mask=vis_div.mask)
        vis_clean_vv=np.ma.array(vis_clean_v,mask=vis_div.mask)
        print np.ma.mean(vis_clean_hh),np.ma.mean(vis_clean_vv)

        d={}
        d['mask']=vis_div.mask
        fs=open(output_file+fname+'_'+str(ant)+'_mask','wb')
        pickle.dump(d,fs,protocol=2)
        fs.close()


        plt.figure(figsize=(14,5.4))
        plt.subplot(121)

        plt.imshow(vis_clean_hh,aspect='auto',extent=(data.freqs[0]/1e6,data.freqs[-1]/1e6,len(timestamps)*2,0))
        plt.xlabel('Freq (MHz)')
        plt.ylabel('time (s)')

        plt.title('clean vis of '+str(fname)+', '+str(ant)+'h', y=1.12)
        plt.colorbar()
        plt.twiny()
        plt.imshow(vis_clean_hh,aspect='auto',extent=(0,len(data.freqs),len(timestamps)*2,0))
        plt.xlabel('channel')

        plt.subplot(122)

        plt.imshow(vis_clean_vv,aspect='auto',extent=(data.freqs[0]/1e6,data.freqs[-1]/1e6,len(timestamps)*2,0))
        plt.xlabel('Freq (MHz)')
        plt.ylabel('time (s)')

        plt.title('clean vis of '+str(fname)+', '+str(ant)+'v', y=1.12)
        plt.colorbar()
        plt.twiny()
        plt.imshow(vis_clean_vv,aspect='auto',extent=(0,len(data.freqs),len(timestamps)*2,0))
        plt.xlabel('channel')

        #plt.savefig(str(fname)+'raw_vis.pdf', bbox_inches='tight')
        plt.savefig(output_file+str(fname)+'_'+ant+'_clean_vis.png', bbox_inches='tight')
        #plt.show()

        #show channel changes

        ch_plot_list=[600,800,1000,2200,2400,2600,2800,3000]
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
                #kv.plot_data(ra[nd_s0],dec[nd_s0], p_data,gsize=60)
                kv.plot_mdata(ra[nd_s0],dec[nd_s0], p_data,gsize=plot_gsize,  grid_method='linear', levels=6, x_mask=1, y_mask=2)
                #plt.plot(p_radec[:,0],p_radec[:,1],'mo')
                if i>len(ch_plot_list)-3:
                    plt.xlabel('R.A.')
                plt.ylabel('Dec')
                plt.title(fname+'_'+recv+'_ch'+str(ch_plot1))
                #plt.plot(p_radec[:,0],p_radec[:,1],'mo')
            plt.savefig(output_file+fname+'_'+recv+'_map.png')
            #plt.show()
    except(Exception):
        print '***skipped becasue of some reason***'
else:
    print '***bad ant, skipped***'

del data
print 'end @ ' + time.asctime(time.localtime(time.time())) +'#'

#END of Level 1 
################################################
#### Jingying Wang <astro.jywang@gmail.com> ####
################################################

