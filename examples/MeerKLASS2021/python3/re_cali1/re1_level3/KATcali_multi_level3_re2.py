#!/usr/bin/env python3
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

print ('start @ ' + time.asctime(time.localtime(time.time())) +'#')

print (katcali.__version__)

print  (plt.rcParams['font.size'], plt.rcParams[u'axes.linewidth'],plt.rcParams['lines.linewidth'])
plt.rcParams['font.size'], plt.rcParams[u'axes.linewidth'],plt.rcParams['lines.linewidth'] = 14, 1.5, 1.5
#plt.rcParams['font.size'], plt.rcParams[u'axes.linewidth'],plt.rcParams['lines.linewidth'] = 10.0, 0.8, 1.5
print  (plt.rcParams['font.size'], plt.rcParams[u'axes.linewidth'],plt.rcParams['lines.linewidth'])
plot_gsize=90

input_file1='/idia/projects/hi_im/raw_vis/MeerKLASS2021/level1/mask/checked/'
input_file2='/idia/projects/hi_im/raw_vis/MeerKLASS2021/level2/data/'
input_file3='../level3/'
output_file='../level3/'
ch_count=0

fname=sys.argv[1]

#select ant, polarization, and one channel to show data calibration
ant=sys.argv[2]
pol=sys.argv[3]
recv=ant+pol
#ch_plot=800

d2=pickle.load(open(input_file2+fname+'_'+str(recv)+'_level2_Tnd_data','rb'),encoding='latin-1')
#print (d2.keys())
assert('Tnd_final' in d2.keys())

data=kio.load_data(fname)
#print data
#print data.obs_script_log

#show the calibrator and bad ants information
target,c0,bad_ants,flux_model=kio.check_ants(fname)


ants_good=[]
for i in np.array(kio.ant_list(data)):
    if i not in bad_ants:
        ants_good.append(i)
    else:
        print (str(i) + ' is bad')
        
print (fname)
print (ants_good)

ch_ref=800
data.select(ants=ant,pol=pol)

corr_id=kio.cal_corr_id(data,recv)
assert(recv==data.corr_products[corr_id][0])
assert(recv==data.corr_products[corr_id][1])
print (corr_id,recv)

vis,flags= kio.call_vis(fname,recv)
vis_backup=vis.copy()
ra,dec,az,el=kio.load_coordinates(data)
timestamps,freqs=kio.load_tf(data)
dump_period=data.dump_period
ang_deg=kio.load_ang_deg(ra,dec,c0)
dp_tt,dp_ss,dp_f,dp_w, dp_t,dp_s,dp_slew,dp_stop=kl.cal_dp_label(data,flags,ant,pol,ch_ref,ang_deg)
dp_sb,dp_se=dp_ss[0],dp_ss[-1]
#dp_u=kl.cal_dp_u(dp_tt,dp_ss)
nd_on_time,nd_cycle,nd_set=kd.cal_nd_basic_para(fname)
print (nd_on_time,nd_cycle,nd_set)
nd_on_edge,nd_off_edge=kd.cal_nd_edges(timestamps,nd_set,nd_cycle,nd_on_time)
print (len(nd_on_edge),len(nd_off_edge))
nd_ratio,nd_0, nd_1x=kd.cal_nd_ratio(timestamps, nd_on_time, nd_on_edge, dump_period)

nd_t0,nd_t1x,nd_s0,nd_s1x,nd_t0_ca,nd_t0_cb,nd_t1x_ca,nd_t1x_cb=kl.cal_label_intersec(dp_tt,dp_ss,nd_0,nd_1x)
labels_1x=kl.cal_label_intersec_complex(dp_tt,dp_ss,nd_0,nd_1x,nd_ratio)

ra1=ra.copy()
for i in range(len(ra)):
    if ra[i]>180:
        ra1[i]=ra[i]-360
ra=ra1

p_radec=np.loadtxt('radio_source2021.txt')

ra_a4059,dec_a4059=-0.74042,-34.76056

p = SkyCoord(data.ra*u.deg,  data.dec*u.deg, frame='icrs')
ang_lim=.5

dp_ptr_list=kl.cal_ptr_mask(p,p_radec,nd_s0, dp_sb,dp_se,ang_lim)


#check with .py result
try:
    d3 = pickle.load(open(input_file1+fname+'_'+ant+'_mask2','rb'))
    print ('# mask2 loaded')
except(Exception):
    d3 = pickle.load(open(input_file1+fname+'_'+ant+'_mask','rb'))
    print ('# mask loaded')
                          
mask_inter=d3['mask']
vis_clean=np.ma.array(vis,mask=mask_inter)

####prepare for data storage#################
gt_param=[None for i in range(len(freqs))]
sm_param=[None for i in range(len(freqs))]
T_map=np.ma.array(np.zeros_like(vis),mask=True)
Tresi_map=np.ma.array(np.zeros_like(vis),mask=True)
Tresi_recover_map=np.ma.array(np.zeros_like(vis),mask=True)
gain_map=np.ma.array(np.zeros_like(vis),mask=True)
Tsm_map=np.ma.array(np.zeros_like(vis),mask=True)
Tel_map=np.ma.array(np.zeros_like(vis),mask=True)

d={}

#Galactic model

nside=64 #healpix nside, 64: Mean Spacing (deg) is 0.9161
#gal_ori=km.cal_Gal_model_np(vis, freqs, ra, dec, ch_plot, ch_plot+1, nside)
gal_ori=km.cal_Gal_model_np2(vis, freqs, ra, dec, 0, len(freqs), nside, model_key=-1)
print ('#Gal model is from Halsam!!!')
gal_ori.flags.writeable=False #avoid change by mistake

#Galactic model from Tsky_block #2023
dg=pickle.load(open(input_file3+fname+'_'+str(ant)+'_Tsky_block_data','rb'))
print (dg.keys())
# assert((dg['ra']==ra).all()) #some diff=360
assert((dg['dec']==dec).all())
Tsky_block=dg['Tsky_block']
gal_ori2=Tsky_block-Tcmb

#spill model 
#Tspill=km.cal_Tspill(el,pol,freqs, ch_plot,2) #fixed version
Tspill_func=km.cal_Tspill_func(el,pol,freqs)

Trec_list=km.cal_Trec(data,ant,pol,freqs)

Tnd_list=d2['Tnd_final']
Tnd_ref_list=d2['Tnd_ref_list']

###raw vis preparsion
#vis_clean_tt=vis_clean.copy()
vis_clean_tt=np.ma.array(vis_clean,mask=gal_ori2.mask) #to add the new count mask, 20230531
print ('# apply the count mask of last round to vis_clean_tt')
vis_clean_tt.mask[:dp_sb,:]=True
vis_clean_tt.mask[dp_se+1:,:]=True

for ch_plot in list(range(550,1051)) + list(range(2150,3101)):
    check_ratio=float(np.array(vis_clean.mask[:,ch_plot]==False).sum())/len(timestamps)
    print (ch_plot, check_ratio)
    if check_ratio>0.3:
        try:
            ###set input params
            Tnd=Tnd_list[ch_plot]
            
            if isinstance(Tnd,np.float)==True:
                print (Tnd)
                #atmosphere emission model
                Tatmo=km.calc_atmosphere_model_1ch(data,ch_plot)
                
                #spill model 
                Tspill=Tspill_func((el,freqs[ch_plot]/1e6))[:,0]

                #elevation related emission model
                Tel=Tspill+Tatmo 
                
                Tgal1=gal_ori[:,ch_plot]
                Tgal2=gal_ori2[:,ch_plot]
                Tgal=Tgal1.copy()
                label_dp=np.ma.where(Tgal2>0)[0]
                print (np.shape(label_dp))
                Tgal[label_dp]=Tgal2[label_dp]
                
                print ('# Tptr has been included in Tgal')

                ####param0
                g0=10.
                Tptr=0 #no point source
                eta_p0=1.0
                Trec0=Trec_list[ch_plot]
                func_sm_param0=[Trec0,0,0,0]
                func_gt_param0=[g0,0,0,0,0]#must be [-5:] from func_obj_sm

                ##fitting
                instru_p=ks.solve_params_sm_v3(timestamps, vis_clean_tt, ch_plot, nd_ratio, Tptr, eta_p0, Tnd, Tel, Tgal,
                                      func_gt_param0, func_sm_param0, nd_0, nd_1x, nd1_weight_plus=0) #modified 2023.4.4

                ###output
                eta_p=instru_p[0]
                sm=instru_p[1:-5]
                gt=instru_p[-5:] 
                print (eta_p, sm, gt)

                m=ks.calc_total_model_sm_v3(timestamps, nd_ratio, Tptr, eta_p, Tnd, Tel, Tgal,  gt, sm, nd_0, nd_1x)
                gain=ks.func_gt(timestamps,gt)
                Tsm=ks.func_sm(timestamps,sm)
                #T=vis_clean[:,ch_plot]/gain
                #residual=(vis_clean[:,ch_plot]-m)/gain
                T=vis_clean_tt[:,ch_plot]/gain #20230604
                residual=(vis_clean_tt[:,ch_plot]-m)/gain
                residual_recover=residual+Tgal-Tgal1
                assert((abs(T[nd_s0]-m[nd_s0]/gain[nd_s0]-residual[nd_s0])<1e-10).all()==True)
                print ('residual in ch'+str(ch_plot)+': '+str(residual[nd_s0].mean())+'+/-'+str(residual[nd_s0].std()))

                ####data need to storage######
                sm_param[ch_plot]=sm
                gt_param[ch_plot]=gt
                T_map[nd_s0,ch_plot]=T[nd_s0]
                Tresi_map[nd_s0,ch_plot]=residual[nd_s0]
                Tresi_recover_map[nd_s0,ch_plot]=residual_recover[nd_s0]
                gain_map[dp_sb:dp_se+1,ch_plot]=gain[dp_sb:dp_se+1]
                Tsm_map[dp_sb:dp_se+1,ch_plot]=Tsm[dp_sb:dp_se+1]
                Tel_map[:,ch_plot]=Tel
                ch_count+=1
                
                if ch_plot==800:

                    #residual map in different resolution
                    plt.figure(figsize=(20,7))
                    plt.subplots_adjust(wspace=0.1,hspace=.4)
                    ax=plt.subplot(221)
                    p_data=T[nd_s0]
                    plt.scatter(ra[nd_s0],dec[nd_s0], c=p_data, vmin=p_data.min(),vmax=p_data.max(), s=8, cmap=kv.cmap1())
                    ax.invert_xaxis()
                    plt.title('calibrated Tmap of '+str(fname)+', '+str(recv)+ ' @'+ str(int(freqs[ch_plot]/1e6))+' MHz')
                    if fname=='1551055211':
                        plt.title('$T_{cal}$ map at '+str(int(freqs[ch_plot]/1e6))+' MHz, '+recv+' of obs190225')
                    plt.xlabel('R.A. (J2000)')
                    plt.ylabel('Dec (J2000)')
                    clb = plt.colorbar()
                    clb.set_label('Kelvin', labelpad=-35, y=1.1, rotation=0)
                    ax=plt.subplot(222)
                    p_data=residual[nd_s0]
                    plt.scatter(ra[nd_s0],dec[nd_s0], c=p_data, vmin=p_data.min(),vmax=p_data.max(), s=8, cmap=kv.cmap1())
                    ax.invert_xaxis()
                    plt.title('residual Tmap of '+str(fname)+', '+str(recv)+ ' @'+ str(int(freqs[ch_plot]/1e6))+' MHz')
                    if fname=='1551055211':
                        plt.title('$T_{res}$ map at '+str(int(freqs[ch_plot]/1e6))+' MHz')
                    plt.xlabel('R.A. (J2000)')
                    plt.ylabel('Dec (J2000)')
                    clb = plt.colorbar()
                    clb.set_label('Kelvin', labelpad=-35, y=1.1, rotation=0)
                    ax=plt.subplot(223)
                    p_data=T[nd_s0]
                    kv.plot_mdata(ra[nd_s0],dec[nd_s0], p_data,gsize=plot_gsize, grid_method='linear', levels=6, x_mask=1, y_mask=2, cmap=kv.cmap1())
                    ax.invert_xaxis()
                    plt.plot(p_radec[:,0],p_radec[:,1],'mo')
                    plt.title('smoothed map of above')
                    #plt.text(146.3, 8.5, 'Kelvin', rotation=0)
                    plt.xlabel('R.A. (J2000)')
                    plt.ylabel('Dec (J2000)')
                    ax=plt.subplot(224)
                    p_data=residual[nd_s0]
                    kv.plot_mdata(ra[nd_s0],dec[nd_s0], p_data,gsize=plot_gsize, grid_method='linear', levels=6, x_mask=1, y_mask=2, cmap=kv.cmap1())
                    ax.invert_xaxis()
                    plt.plot(p_radec[:,0],p_radec[:,1],'mo')
                    plt.title('smoothed map of above')
                    #plt.text(146.3, 8.5, 'Kelvin', rotation=0)
                    plt.xlabel('R.A. (J2000)')
                    plt.ylabel('Dec (J2000)')
                    #plt.savefig('F_calibrated_map.pdf', bbox_inches='tight')
                    plt.savefig(output_file+'F_'+fname+'_'+recv+'_calibrated_map.png', bbox_inches='tight')
                    #plt.show()

            else:              
                print ("* no Tnd value for level-3 calibration")
        except(Exception):
            print ('***channel'+ str(ch_plot) +' failed')
    else:
        print ('***channel'+ str(ch_plot) +' skipped')
####save data####
d['gt_param']=gt_param
d['sm_param']=sm_param
d['T_map']=T_map
d['Tresi_map']=Tresi_map
d['Tresi_recover_map']=Tresi_recover_map
d['gain_map']=gain_map
d['Tel_map']=Tel_map
d['Tgal_map']=gal_ori
d['Tsm_map']=Tsm_map
d['nd_s0']=nd_s0
d['timestamps']=timestamps
d['ra']=ra
d['dec']=dec
fs=open(output_file+fname+'_'+recv+'_level3_re1_data','wb')
pickle.dump(d,fs,protocol=2)
fs.close()

print ('# total fitted channel number: '+str(ch_count))
print ('end @ ' + time.asctime(time.localtime(time.time())) +'#')
print ('Hakuna Matata')
################################################
#### Jingying Wang <astro.jywang@gmail.com> ####
################################################
