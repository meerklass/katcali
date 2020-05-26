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

output_file='level2_output/'

fname=sys.argv[1]
data=kio.load_data(fname)
#print data

#show the calibrator and bad ants information
target,c0,bad_ants,flux_model=kio.check_ants(fname)
# Select ant and polarization, then load data in 
ant=sys.argv[2]
pol=sys.argv[3]

print '***Level2 calibartion for '+fname+' '+ant+pol+' start...'

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
T_map=np.ma.array(np.zeros_like(vis),mask=True)
Tresi_map=np.ma.array(np.zeros_like(vis),mask=True)
gain_map=np.ma.array(np.zeros_like(vis),mask=True)
Tel_map=np.ma.array(np.zeros_like(vis),mask=True)

Tnd_ref_list=[None for i in range(len(freqs))]
Tnda_list=[None for i in range(len(freqs))]
Tndb_list=[None for i in range(len(freqs))]

gta_param=[None for i in range(len(freqs))]
gtb_param=[None for i in range(len(freqs))]

sma_param=[None for i in range(len(freqs))]
smb_param=[None for i in range(len(freqs))]

ratioa_list=[None for i in range(len(freqs))]
ratiob_list=[None for i in range(len(freqs))]
d={}

mask_nd_t0=mask_inter.copy()
mask_nd_t0[:,:]=True
assert((mask_nd_t0==True).all()==True)

####prepare data for multi-band calibration####
#beam
beam_select='me'
Npix=513
Ddeg=5
pattern_fband=kb.load_pattern_fband(beam_select,pol)
Aeff_max_fband=kb.load_Aeff_max_fband(beam_select,pol)
x_pix,y_pix=kb.cal_pix_params(data,c0,Npix,Ddeg)
#Galactic model
nside=64 #healpix nside, 64: Mean Spacing (deg) is 0.9161
gal_ori=km.cal_Gal_model_np(vis, freqs, ra, dec, 0, len(freqs), nside)
gal_ori.flags.writeable=False #avoid change by mistake
gal=gal_ori.copy() #will change for some track data
#Spill
Tspill_func=km.cal_Tspill_func(el,pol,freqs)
#diode noise
Tnd_std,Tnd_spl=km.Tnd_spl(data, ant,pol)
print Tnd_std
#receiver
Trec_list=km.cal_Trec(data,ant,pol,freqs)

#########################################################################
ch_list0=range(550,1051) + range(2150,3101) ##set channels
ch_list=ch_list0[:]

####check data exists
try:
    d0=pickle.load(open(output_file+str(fname)+'_'+str(recv)+'_level2_data_tmp'))
    T_map=d0['T_map']
    Tresi_map=d0['Tresi_map']
    gain_map=d0['gain_map']
    Tel_map=d0['Tel_map']
    Tnd_ref_list=d0['Tnd_ref_list']
    Tnda_list=d0['Tnda_list']
    Tndb_list=d0['Tndb_list']
    gta_param=d0['gta_param']
    gtb_param=d0['gtb_param']
    sma_param=d0['sma_param']
    smb_param=d0['smb_param']
    ratioa_list=d0['ratioa_list']
    ratiob_list=d0['ratiob_list']
    mask_nd_t0=d0['mask_nd_t0']
    print '#### pre data loaded!!!'
    for ch_i in ch_list0:
        if (mask_nd_t0[:,ch_i]==True).all()==False:
            ch_list.remove(ch_i)
            print 'ch'+str(ch_i)+' data exists.'
except IOError:
    print '#### no pre data exists.'

#########################################################################
#select ant, polarization, and one channel to show data calibration
for ch_plot in ch_list:  #range(550,1051) + range(2150,3101):##set channels
    check_ratio=float(np.array(vis_clean.mask[:,ch_plot]==False).sum())/len(timestamps)
    print ch_plot, check_ratio
    if check_ratio>0.6:
        try:
            assert(np.shape(data)[2]==1)
            dp_tt,dp_ss,dp_f,dp_t,dp_s=kl.label_dump_1ch(data,ant,pol,flags,ch_plot)
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
            #BM-III:pattern
            if fname in ['1551055211','1551037708', '1579725085', '1580260015']:
                dp_ca,dp_cb,dp_c0a, dp_c1a,dp_c2a,dp_c3a,dp_c4a,dp_c0b,dp_c1b,dp_c2b,dp_c3b,dp_c4b=kl.cal_dp_c(fname,data,ant,pol,flags,ch_plot,dp_tt,dp_ss,ang_deg)
            else:
                dp_ca,dp_cb,dp_c0a, dp_c1a,dp_c0b,dp_c1b=kl.cal_dp_c(fname,data,ant,pol,flags,ch_plot,dp_tt,dp_ss,ang_deg)
            #calculate position
            T_ptr2,pattern,x_pix_max,y_pix_max=kb.cal_BMIII_1ch(data,ch_plot,flux_model, dp_ca,dp_cb,pattern_fband,x_pix,y_pix,Aeff_max_fband)
            #####choose beam model
            T_ptr=T_ptr2 #BM-III 
            #Galactic model
            if fname in ['1551037708', '1551055211']:
                dp_gal=list(dp_c0a)+list(dp_c0b)
            else:
                dp_gal=list(dp_c1a)+list(dp_c1b)

            gal_reset=np.mean(gal_ori[dp_gal,ch_plot])
            gal[dp_tt,ch_plot]=np.ones_like(gal_ori[dp_tt,ch_plot])*gal_reset
            Tgal=gal[:,ch_plot]

            ### track before scan#
            #select raw vis for track befor scan
            visa_ptr = vis_clean.copy()
            for i in range(len(timestamps)):
                if i not in dp_ca:
                    visa_ptr.mask[i,:]=True
            ####set input parameters
            ga0,gb0=ks.cal_gain0(fname,data,ant,pol,flags,ch_plot,dp_tt,dp_ss,ang_deg,T_ptr,vis_clean) #gain level
            print ga0, gb0
            assert(isinstance(ga0,np.float))
            assert(isinstance(gb0,np.float))
            
            Trec0=Trec_list[ch_plot]
            eta_p0=1.0
            func_sm_param0=[Trec0]
            func_gt_param0=[ga0,0,0,0,0]
            ratio0=0.5
            ####fitting
            instru_pa=ks.solve_params0(timestamps, visa_ptr, ch_plot, nd_ratio, ratio0, T_ptr, eta_p0, Tnd_ref, Tnd_std, Tel, Tgal,
                                  func_gt_param0, func_sm_param0, nd_0, nd_1a, nd_1b)
            ####get fitting result
            Tnda=instru_pa[0]
            eta_pa=instru_pa[1]
            sma=instru_pa[2]
            gta=instru_pa[3:-1]
            ratioa=instru_pa[-1]
            print Tnda, eta_pa, ratioa , sma, gta

            ### track after scan#
            #select raw vis for track befor scan
            visb_ptr = vis_clean.copy()
            for i in range(len(timestamps)):
                if i not in dp_cb:
                    visb_ptr.mask[i,:]=True
            ####set input parameters
            eta_p0=1.0
            func_sm_param0=[Trec0]
            func_gt_param0=[gb0,0,0,0,0]
            ratio0=0.5
            ####fitting######
            instru_pb=ks.solve_params0(timestamps, visb_ptr, ch_plot, nd_ratio, ratio0, T_ptr, eta_p0, Tnd_ref, Tnd_std, Tel, Tgal,
                                  func_gt_param0, func_sm_param0, nd_0, nd_1a, nd_1b)
            ######get fitting result#####
            Tndb=instru_pb[0]
            eta_pb=instru_pb[1]
            smb=instru_pb[2]
            gtb=instru_pb[3:-1]
            ratiob=instru_pb[-1]
            print Tndb, eta_pb, ratiob , smb, gtb

            #result
            m=ks.calc_total_model(timestamps, nd_ratio, ratioa, T_ptr, eta_pa, Tnda, Tel, Tgal, gta, sma, nd_0, nd_1a, nd_1b)
            ma=np.ma.array(m,mask=visa_ptr[:,ch_plot].mask)
            ga=ks.func_gt(timestamps,gta)
            resia=(visa_ptr[:,ch_plot]-ma)/ga

            m=ks.calc_total_model(timestamps, nd_ratio, ratiob, T_ptr, eta_pb, Tndb, Tel, Tgal, gtb, smb, nd_t0, nd_t1a, nd_t1b)
            mb=np.ma.array(m,mask=visb_ptr[:,ch_plot].mask)
            gb=ks.func_gt(timestamps,gtb)
            resib=(visb_ptr[:,ch_plot]-mb)/gb

            ####data need to storage######
            Tnd_ref_list[ch_plot]=Tnd_ref
            Tnda_list[ch_plot]=Tnda
            Tndb_list[ch_plot]=Tndb

            gta_param[ch_plot]=gta
            gtb_param[ch_plot]=gtb

            gain_map[dp_ca,ch_plot]=ga[dp_ca]
            gain_map[dp_cb,ch_plot]=gb[dp_cb]

            calT_tra=visa_ptr[:,ch_plot]/ga
            calT_trb=visb_ptr[:,ch_plot]/gb

            assert((abs(calT_tra[dp_ca]-ma[dp_ca]/ga[dp_ca]-resia[dp_ca])<1e-10).all()==True)
            assert((abs(calT_trb[dp_cb]-mb[dp_cb]/gb[dp_cb]-resib[dp_cb])<1e-10).all()==True)
            T_map[dp_ca,ch_plot]=calT_tra[dp_ca]
            T_map[dp_cb,ch_plot]=calT_trb[dp_cb]

            Tresi_map[dp_ca,ch_plot]=resia[dp_ca]
            Tresi_map[dp_cb,ch_plot]=resib[dp_cb]

            sma_param[ch_plot]=sma
            smb_param[ch_plot]=smb

            ratioa_list[ch_plot]=ratioa
            ratiob_list[ch_plot]=ratiob

            Tel_map[:,ch_plot]=Tel

            mask_nd_t0[nd_t0,ch_plot]=False

            if ch_plot % 25 == 0:
                d['T_map']=T_map
                d['Tresi_map']=Tresi_map
                d['gain_map']=gain_map
                d['Tel_map']=Tel_map
                d['Tnd_ref_list']=Tnd_ref_list
                d['Tnda_list']=Tnda_list
                d['Tndb_list']=Tndb_list
                d['gta_param']=gta_param
                d['gtb_param']=gtb_param
                d['sma_param']=sma_param
                d['smb_param']=smb_param
                d['ratioa_list']=ratioa_list
                d['ratiob_list']=ratiob_list
                d['mask_nd_t0']=mask_nd_t0
                d['timestamps']=timestamps
                d['ra']=ra
                d['dec']=dec
                fs=open(output_file+str(fname)+'_'+str(recv)+'_level2_data_tmp','wb')
                pickle.dump(d,fs,protocol=2)
                fs.close()
                print '*data saved at ch'+str(ch_plot)
        except(Exception):
            print '***channel'+ str(ch_plot) +' failed'
    else:
        print '***channel'+ str(ch_plot) +' skipped'        

####save data####
d['T_map']=T_map
d['Tresi_map']=Tresi_map
d['gain_map']=gain_map
d['Tel_map']=Tel_map
d['Tnd_ref_list']=Tnd_ref_list
d['Tnda_list']=Tnda_list
d['Tndb_list']=Tndb_list
d['gta_param']=gta_param
d['gtb_param']=gtb_param
d['sma_param']=sma_param
d['smb_param']=smb_param
d['ratioa_list']=ratioa_list
d['ratiob_list']=ratiob_list
d['mask_nd_t0']=mask_nd_t0
d['timestamps']=timestamps
d['ra']=ra
d['dec']=dec
fs=open(output_file+str(fname)+'_'+str(recv)+'_level2_data','wb')
pickle.dump(d,fs,protocol=2)
fs.close()


d2={}
d2['Tnd_ref_list']=Tnd_ref_list
d2['Tnda_list']=Tnda_list
d2['Tndb_list']=Tndb_list
fs=open(output_file+str(fname)+'_'+str(recv)+'_level2_Tnd_data','wb')
pickle.dump(d2,fs,protocol=2)
fs.close()

print '***Level2 calibartion for '+fname+' '+ant+pol
print 'end @ ' + time.asctime(time.localtime(time.time())) +'#'
print 'Hakuna Matata' 
################################################
#### Jingying Wang <astro.jywang@gmail.com> ####
################################################
