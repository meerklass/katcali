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


print 'start @ ' + time.asctime(time.localtime(time.time())) +'#'
print katcali.__version__

output_file='level3_output/'

fname=sys.argv[1]

data=kio.load_data(fname)
#print data

#show the calibrator and bad ants information
target,c0,bad_ants,flux_model=kio.check_ants(fname)
# Select ant and polarization, then load data in 
ant=sys.argv[2]
pol=sys.argv[3]

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
ang_deg=kio.load_ang_deg(ra,dec,c0)
nd_1a,nd_1b,nd_1,nd_0=kd.call_nd_1_list(fname,timestamps)

# RFI flagging
try:
    d3 = pickle.load(open('/idia/projects/hi_im/raw_vis/katcali_output/level1_output/mask/'+fname+'_'+ant+'_mask2'))
    print 'mask2 loaded'
except(Exception):
    d3 = pickle.load(open('/idia/projects/hi_im/raw_vis/katcali_output/level1_output/mask/'+fname+'_'+ant+'_mask'))
    print 'mask loaded'
mask_inter=d3['mask']
vis_clean=np.ma.array(vis,mask=mask_inter)

####prepare for data storage#################
gt_param=[None for i in range(len(freqs))]
sm_param=[None for i in range(len(freqs))]
ratio_list=[None for i in range(len(freqs))]

T_map=np.ma.array(np.zeros_like(vis),mask=True)
Tresi_map=np.ma.array(np.zeros_like(vis),mask=True)
gain_map=np.ma.array(np.zeros_like(vis),mask=True)
Tsm_map=np.ma.array(np.zeros_like(vis),mask=True)
Tel_map=np.ma.array(np.zeros_like(vis),mask=True)
Tgal_map=np.ma.array(np.zeros_like(vis),mask=True)

mask_nd_s0=mask_inter.copy()
mask_nd_s0[:,:]=True
assert((mask_nd_s0==True).all()==True)

d={}
####prepare data for multi-band calibration####
#Galactic model
nside=64 #healpix nside, 64: Mean Spacing (deg) is 0.9161
gal_ori=km.cal_Gal_model_np(vis, freqs, ra, dec, 0, len(freqs), nside)
gal_ori.flags.writeable=False #avoid change by mistake
#Spill
Tspill_func=km.cal_Tspill_func(el,pol,freqs)
#diode noise
Tnd_std,Tnd_spl=km.Tnd_spl(data, ant,pol)
print Tnd_std
#receiver
Trec_list=km.cal_Trec(data,ant,pol,freqs)
#Tnd
d1=pickle.load(open('/scratch/users/jywang/WiggleZ_11hr/level2_output/done/'+str(fname)+'_'+str(recv)+'_level2_Tnd_data'))
#print d.keys()
Tnda_list=d1['Tnda_list']
Tndb_list=d1['Tndb_list']

#select ant, polarization, and one channel to show data calibration
for ch_plot in range(550,1051) + range(2150,3101):##set channels

    check_ratio=float(np.array(vis_clean.mask[:,ch_plot]==False).sum())/len(timestamps)
    print ch_plot, check_ratio
    if check_ratio>0.6:
        try:

            dp_tt,dp_ss,dp_f,dp_t,dp_s=kl.label_dump_1ch(data,ant,pol,flags,ch_plot)
            assert(np.shape(data)[2]==1)
            dp_sb=dp_ss[0]
            dp_se=dp_ss[-1]
            nd_s1a,nd_s1b,nd_s1,nd_s0=kd.cal_nds_list(dp_ss,nd_1a,nd_1b,nd_1,nd_0)#dp_ss here, not dp_s
            nd_t1a,nd_t1b,nd_t1,nd_t0=kd.cal_ndt_list(dp_tt,nd_1a,nd_1b,nd_1,nd_0)#dp_tt here, not dp_t

            ## load the foreground models
            #spill model 
            Tspill=Tspill_func((el,freqs[ch_plot]/1e6))[:,0]
            Tatmo=km.calc_atmosphere_model_1ch(data,ch_plot)
            Tel=Tspill+Tatmo 
            #load the diode injection model and get a reference value
            Tnd_ref=Tnd_spl(freqs[ch_plot]/1e9)

            #Galactic model
            Tgal=gal_ori[:,ch_plot] #different to track part

            #Tnd
            Tnda=Tnda_list[ch_plot]
            Tndb=Tndb_list[ch_plot]
            if fname in ['1551037708']:
                Tnd=Tndb
            else:
                Tnd=(Tnda+Tndb)/2.
            print Tnda,Tndb,Tnd
            assert(isinstance(Tnd,np.float))
            
            ###raw vis preparsion
            vis_clean_ss=vis_clean.copy()
            vis_clean_ss.mask[:dp_sb,:]=True
            vis_clean_ss.mask[dp_se+1:,:]=True

            ####param0
            g0=10.
            Tptr=0 #no point source
            eta_p0=1.0
            Trec0=Trec_list[ch_plot]
            func_sm_param0=[Trec0,0,0,0]
            func_gt_param0=[g0,0,0,0,0]#must be [-6:-1] from func_obj_sm
            ratio0=0.5

            ##fitting
            instru_p=ks.solve_params_sm(timestamps, vis_clean_ss, ch_plot, nd_ratio, ratio0, Tptr, eta_p0, Tnd, Tel, Tgal,
                                  func_gt_param0, func_sm_param0, nd_0, nd_1a, nd_1b)

            ###output
            eta_p=instru_p[0]
            sm=instru_p[1:-6]
            Tsm=ks.func_sm(timestamps,sm)
            gt=instru_p[-6:-1] #must be [-5:-1] from func_obj_sm
            gain=ks.func_gt(timestamps,gt)
            ratio=instru_p[-1]
            print eta_p, ratio, sm, gt

            T=vis_clean[:,ch_plot]/gain
            m=ks.calc_total_model_sm(timestamps, nd_ratio, ratio, Tptr, eta_p, Tnd, Tel, Tgal,  gt, sm, nd_0, nd_1a, nd_1b)
            residual=(vis_clean[:,ch_plot]-m)/gain
            assert((abs(T[nd_s0]-m[nd_s0]/gain[nd_s0]-residual[nd_s0])<1e-10).all()==True)

            #print gain[dp_sb:dp_se+1].mean(),gain[dp_sb:dp_se+1].std()
            #print m[dp_sb:dp_se+1].mean(),m[dp_sb:dp_se+1].std()
            #print T[nd_s0].mean(),T[nd_s0].std()
            #print residual[nd_s0].mean(),residual[nd_s0].std()

            ####data need to storage######

            sm_param[ch_plot]=sm
            gt_param[ch_plot]=gt
            ratio_list[ch_plot]=ratio

            T_map[nd_s0,ch_plot]=T[nd_s0]
            Tresi_map[nd_s0,ch_plot]=residual[nd_s0]
            gain_map[dp_sb:dp_se+1,ch_plot]=gain[dp_sb:dp_se+1]
            Tsm_map[dp_sb:dp_se+1,ch_plot]=Tsm[dp_sb:dp_se+1]
            Tel_map[:,ch_plot]=Tel

            mask_nd_s0[nd_s0,ch_plot]=False

            if ch_plot % 50 == 0:

                d['gt_param']=gt_param
                d['sm_param']=sm_param
                d['ratio_list']=ratio_list

                d['T_map']=T_map
                d['Tresi_map']=Tresi_map
                d['gain_map']=gain_map
                d['Tel_map']=Tel_map
                #d['Tgal_map']=Tgal_map
                d['Tsm_map']=Tsm_map

                d['mask_nd_s0']=mask_nd_s0

                d['timestamps']=timestamps
                d['ra']=ra
                d['dec']=dec
                fs=open(output_file+str(fname)+'_'+str(recv)+'_level3_data_tmp','wb')
                pickle.dump(d,fs,protocol=2)
                fs.close()

           
        except(Exception):
            print '***channel'+ str(ch_plot) +' failed'
    else:
        print '***channel'+ str(ch_plot) +' skipped'        
          
####save data####
Tgal_map=gal_ori.copy()

d['gt_param']=gt_param
d['sm_param']=sm_param
d['ratio_list']=ratio_list

d['T_map']=T_map
d['Tresi_map']=Tresi_map
d['gain_map']=gain_map
d['Tel_map']=Tel_map
d['Tgal_map']=Tgal_map
d['Tsm_map']=Tsm_map

d['mask_nd_s0']=mask_nd_s0

d['timestamps']=timestamps
d['ra']=ra
d['dec']=dec
fs=open(output_file+str(fname)+'_'+str(recv)+'_level3_data','wb')
pickle.dump(d,fs,protocol=2)
fs.close()

print 'end @ ' + time.asctime(time.localtime(time.time())) +'#'
print 'Hakuna Matata' 
################################################
#### Jingying Wang <astro.jywang@gmail.com> ####
################################################
