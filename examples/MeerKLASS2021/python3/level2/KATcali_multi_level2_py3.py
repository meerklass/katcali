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
from matplotlib.colors import LogNorm
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

plt.rcParams['font.size'], plt.rcParams[u'axes.linewidth'],plt.rcParams['lines.linewidth'] = 14, 1.5, 1.5

input_file='/idia/projects/hi_im/raw_vis/MeerKLASS2021/level1/mask/checked/'
#output_file='./'
output_file='/scratch3/users/jywang/MeerKLASS2021/level2/'

d={}
d2={}
count_ch=0

fname=sys.argv[1]
ant=sys.argv[2]
pol=sys.argv[3]

data=kio.load_data(fname)
#print data
target,c0,bad_ants,flux_model=kio.check_ants(fname)
if target=='PKS1934-638':
    tag='1934-638'
if target=='PictorA':
    tag='PictorA'

ants_good=[]
for i in np.array(kio.ant_list(data)):
    if i not in bad_ants:
        ants_good.append(i)
    else:
        print (str(i) + ' is bad')
        
print (fname)
print (ants_good)

if ant in ants_good:
    
    data.select()
    data.select(targets=tag+'_u0.5')
    target_start=data.target_indices[0]
    data.select()


    nd_on_time,nd_cycle,nd_set=kd.cal_nd_basic_para(fname)
    print (nd_on_time,nd_cycle,nd_set)

    # Select ant and polarization, then load data in 
    #select ant, polarization, and one channel to show data calibration

    #load data, labels, and parameters
    ch_ref=800
    data.select(ants=ant,pol=pol)
    recv=ant+pol
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

    az_corr=az.copy()
    for i in range(len(az)):
        if az[i]>180:
            az_corr[i]=az[i]-360

    nd_on_edge,nd_off_edge=kd.cal_nd_edges(timestamps,nd_set,nd_cycle,nd_on_time)
    print (len(nd_on_edge),len(nd_off_edge))
    nd_ratio,nd_0, nd_1x=kd.cal_nd_ratio(timestamps, nd_on_time, nd_on_edge, dump_period)

    # RFI flagging
    try:
        d3 = pickle.load(open(input_file+fname+'_'+ant+'_mask2', 'rb'))
        print ('mask2 loaded')
    except(Exception):
        d3 = pickle.load(open(input_file+fname+'_'+ant+'_mask', 'rb'))
        print ('mask loaded')

    mask_inter=d3['mask']
    #vis_clean=np.ma.array(vis,mask=mask_inter)
    vis_clean=np.ma.array(vis,mask=flags)
    vis_clean=np.ma.array(vis_clean,mask=mask_inter)


    # calibrate the diode noise using point source calibrator 
    ## load the foreground models
    #Galactic model
    nside=64 #healpix nside, 64: Mean Spacing (deg) is 0.9161
    #gal_ori=km.cal_Gal_model_np(vis, freqs, ra, dec, ch_plot, ch_plot+1, nside)
    gal_ori=km.cal_Gal_model_np2(vis, freqs, ra, dec, 0, len(freqs), nside, model_key=-1)
    print ('#Gal model is from Halsam!!!')
    gal_ori.flags.writeable=False #avoid change by mistake
    gal=gal_ori.copy() #will change for some track data

    #select beam pattern model
    beam_select='me'
    ###20221220 new bc BMIII -> BMIII_1ch
    Npix=513
    Ddeg=5
    pattern_fband=kb.load_pattern_fband(beam_select,pol)
    Aeff_max_fband=kb.load_Aeff_max_fband(beam_select,pol)
    x_pix,y_pix=kb.cal_pix_params(data,c0,Npix,Ddeg)

    #BM-III:pattern
    dp_ca,dp_cb,dp_c0a, dp_c1a,dp_c2a,dp_c3a,dp_c4a,dp_c0b,dp_c1b,dp_c2b,dp_c3b,dp_c4b,dp_c1a1,dp_c1a2,dp_c1a3,dp_c1b1,dp_c1b2,dp_c1b3=kl.cal_dp_c(fname,data,ant,pol,flags,ch_ref,dp_tt,dp_ss,ang_deg, target_start=target_start,n_src_off=4)

    #Spill
    Tspill_func=km.cal_Tspill_func(el,pol,freqs)

    #diode noise
    Tnd_std,Tnd_spl=km.Tnd_spl(data, ant,pol)
    Trec_list=km.cal_Trec(data,ant,pol,freqs)

    #load the scan and track labels 
    dp_u=kl.cal_dp_u(dp_tt,dp_ss)

    visa_ptr = vis_clean.copy()
    visb_ptr = vis_clean.copy()
    for i in range(len(timestamps)):
        if i not in dp_ca:
            visa_ptr.mask[i,:]=True
        if i not in dp_cb:
            visb_ptr.mask[i,:]=True
    ####prepare for data storage#################
    Tnd_ref_list=[None for i in range(len(freqs))]
    Tel_map=np.ma.array(np.zeros_like(vis),mask=True)
    
    T_map=np.ma.array(np.zeros_like(vis),mask=True)
    Tresi_map=np.ma.array(np.zeros_like(vis),mask=True)
    gain_map=np.ma.array(np.zeros_like(vis),mask=True)
    Tnda_list=[None for i in range(len(freqs))]
    Tndb_list=[None for i in range(len(freqs))]
    NRMSE1_list=[None for i in range(len(freqs))]
    NRMSE2_list=[None for i in range(len(freqs))]
    gta_param=[None for i in range(len(freqs))]
    gtb_param=[None for i in range(len(freqs))]
    sma_param=[None for i in range(len(freqs))]
    smb_param=[None for i in range(len(freqs))]
    
    if target=='PKS1934-638':
        T_map2=np.ma.array(np.zeros_like(vis),mask=True)
        Tresi_map2=np.ma.array(np.zeros_like(vis),mask=True)
        gain_map2=np.ma.array(np.zeros_like(vis),mask=True)
        Tnda_list2=[None for i in range(len(freqs))]
        Tndb_list2=[None for i in range(len(freqs))]
        NRMSE11_list=[None for i in range(len(freqs))]
        NRMSE22_list=[None for i in range(len(freqs))]
        gta_param2=[None for i in range(len(freqs))]
        gtb_param2=[None for i in range(len(freqs))]
        sma_param2=[None for i in range(len(freqs))]
        smb_param2=[None for i in range(len(freqs))]       
    


    for ch_plot in list(range(550,1051)) + list(range(2150,3101)):

        check_ratio=float(np.array(vis_clean.mask[:,ch_plot]==False).sum())/len(timestamps)
        print (ch_plot, check_ratio)
        if check_ratio>0.3: #0.6
            try:

                #spill model 
                Tspill=Tspill_func((el,freqs[ch_plot]/1e6))[:,0]

                #atmosphere emission model
                Tatmo=km.calc_atmosphere_model_1ch(data,ch_plot)

                #atmosphere emission model
                Tatmo_abs_1ch=km.calc_atmosphere_trans_factor_1ch(data,ch_plot) #old name:calc_atmosphere_abs_factor_1ch

                #elevation related emission model
                Tel=Tspill+Tatmo 

                #load the diode injection model and get a reference value
                #note: diode version are different dish by dish!
                Tnd_ref=Tnd_spl(freqs[ch_plot]/1e9)
                print (Tnd_std,Tnd_ref)

                ## load calibrator model: related to the beam model##################



                #calculate position
                #T_ptr2,pattern,pix_label=kb.cal_BMIII(fname,data,ch_plot,ant,pol,flux_model,c0,dp_ca,dp_cb,ang_deg,beam_select)
                #x_pix,y_pix,x_pix_max,y_pix_max=pix_label
                T_ptr2,pattern,x_pix_max,y_pix_max=kb.cal_BMIII_1ch(data,ch_plot,flux_model, dp_ca,dp_cb,pattern_fband,x_pix,y_pix,Aeff_max_fband)

                Tgal=gal[:,ch_plot]

                ## calibrate diode noise

                #####choose beam model
                T_ptr=T_ptr2 #BM-III 

                ### track before scan########################

                ####set input parameters
                ga0,gb0=ks.cal_gain0(fname,data,ant,pol,flags,ch_plot,dp_tt,dp_ss,ang_deg,T_ptr,vis_clean,n_src_off=4,target_start=target_start) #gain level ##need to check the label still works for 2021 data!!!!!!
                print (ga0,gb0)
                assert(isinstance(ga0,np.float))
                assert(isinstance(gb0,np.float))

                Trec0=Trec_list[ch_plot]
                print (Trec0)
                eta_p0=1.0
                func_sm_param0=[Trec0]
                func_gt_param0=[ga0,0,0,0,0]
                print (Tnd_ref)
                ####fitting
                instru_pa=ks.solve_params0_v3(timestamps, visa_ptr, ch_plot, nd_ratio, T_ptr, eta_p0, Tnd_ref, Tnd_std, Tel, Tgal, func_gt_param0, func_sm_param0, nd_0, nd_1x)

                ####get fitting result
                Tnda=instru_pa[0]
                eta_pa=instru_pa[1]
                sma=instru_pa[2]
                gta=instru_pa[3:]

                print (Tnda, eta_pa, sma, gta)

                ### track after scan########################

                ####set input parameters
                #Trec0=Trec_list[ch_plot]
                print (Trec0)
                eta_p0=1.0
                func_sm_param0=[Trec0]
                func_gt_param0=[gb0,0,0,0,0]

                print (Tnd_ref)

                ####fitting######
                instru_pb=ks.solve_params0_v3(timestamps, visb_ptr, ch_plot, nd_ratio, T_ptr, eta_p0, Tnd_ref, Tnd_std, Tel, Tgal,
                                      func_gt_param0, func_sm_param0, nd_0, nd_1x)

                ######get fitting result#####
                Tndb=instru_pb[0]
                eta_pb=instru_pb[1]
                smb=instru_pb[2]
                gtb=instru_pb[3:]

                print (Tndb, eta_pb, smb, gtb)

                m=ks.calc_total_model_v3(timestamps, nd_ratio, T_ptr, eta_pa, Tnda, Tel, Tgal, gta, sma, nd_0, nd_1x)
                ma=np.ma.array(m,mask=visa_ptr[:,ch_plot].mask)
                m=ks.calc_total_model_v3(timestamps, nd_ratio, T_ptr, eta_pb, Tndb, Tel, Tgal, gtb, smb, nd_0, nd_1x)
                mb=np.ma.array(m,mask=visb_ptr[:,ch_plot].mask)
                ga=ks.func_gt(timestamps,gta)
                resia=(visa_ptr[:,ch_plot]-ma)/ga
                gb=ks.func_gt(timestamps,gtb)
                resib=(visb_ptr[:,ch_plot]-mb)/gb

                NRMSE1=ks.cal_NRMSE(ma/ga,resia)
                NRMSE2=ks.cal_NRMSE(mb/gb,resib)

                if ch_plot==800:
                    ##show model and raw vis
                    plt.figure(figsize=(20,8))
                    plt.subplots_adjust(wspace=0.1,hspace=0.25)
                    plt.subplot(221)

                    plt.plot(timestamps[dp_c0a]-timestamps[0],ma[dp_c0a],'r-',drawstyle='steps-mid')
                    plt.plot(timestamps[dp_c0a]-timestamps[0],visa_ptr[dp_c0a,ch_plot],'g--',drawstyle='steps-mid')
                    plt.plot(timestamps[dp_c1a]-timestamps[0],ma[dp_c1a],'r-',drawstyle='steps-mid')
                    plt.plot(timestamps[dp_c1a]-timestamps[0],visa_ptr[dp_c1a,ch_plot],'g--',drawstyle='steps-mid')
                    #if fname in ['1551055211','1551037708','1630519596']:
                    plt.plot(timestamps[dp_c2a]-timestamps[0],ma[dp_c2a],'r-',drawstyle='steps-mid')
                    plt.plot(timestamps[dp_c3a]-timestamps[0],ma[dp_c3a],'r-',drawstyle='steps-mid')
                    plt.plot(timestamps[dp_c4a]-timestamps[0],ma[dp_c4a],'r-',drawstyle='steps-mid')
                    plt.plot(timestamps[dp_c2a]-timestamps[0],visa_ptr[dp_c2a,ch_plot],'g--',drawstyle='steps-mid')
                    plt.plot(timestamps[dp_c3a]-timestamps[0],visa_ptr[dp_c3a,ch_plot],'g--',drawstyle='steps-mid')
                    plt.plot(timestamps[dp_c4a]-timestamps[0],visa_ptr[dp_c4a,ch_plot],'g--',drawstyle='steps-mid')
                    #plt.xlabel('time (s)')
                    plt.ylabel('raw signal')
                    plt.legend(['model','raw signal'],ncol=2)
                    plt.ylim(np.nanmin(visa_ptr[dp_c0a,ch_plot])-20,np.nanmax(visa_ptr[dp_c1a,ch_plot])+30)
                    plt.title('Fitting curve of track-I at '+ str(int(freqs[ch_plot]/1e6)) +' MHz, '+str(fname)+', '+str(recv))
                    if fname=='1551055211':
                        plt.title('Fitting curve of track-I at '+ str(int(freqs[ch_plot]/1e6)) +' MHz, '+str(recv)+' of obs190225')
                    plt.subplot(223)

                    plt.plot([timestamps[dp_c0a][0]-timestamps[0],timestamps[dp_c0a][-1]-timestamps[0]],[0,0],'k-.')
                    plt.plot(timestamps[dp_c0a]-timestamps[0],resia[dp_c0a],'.-')
                    plt.plot(timestamps[dp_c1a]-timestamps[0],resia[dp_c1a],'.-')
                    #if fname in ['1551055211','1551037708','1630519596']:
                    plt.plot(timestamps[dp_c2a]-timestamps[0],resia[dp_c2a],'.-')
                    plt.plot(timestamps[dp_c3a]-timestamps[0],resia[dp_c3a],'.-')
                    plt.plot(timestamps[dp_c4a]-timestamps[0],resia[dp_c4a],'.-')
                    plt.plot([timestamps[dp_c0a][0]-timestamps[0],timestamps[dp_c4a][-1]-timestamps[0]],[0,0],'k-.')
                    plt.xlabel('time (s)')
                    plt.ylabel('Tres (K)')
                    plt.ylim(-0.15,0.15)
                    #plt.legend(['line of T=0 K','track calibrator outskirt', 'track calibrator center'],ncol=3)
                    #if fname in ['1551055211','1551037708','1579725085', '1580260015','1630519596']:
                    plt.legend(['zero line','track c0', 'track c1', 'track c2', 'track c3','track c4'], fontsize=12, ncol=3)
                    plt.ylim(-0.19,0.19)
                    plt.title('$T_{res}$ of the fitted curve above')
                    #plt.savefig('caliA_ch'+str(ch_plot)+'.png')
                    ##show model and raw vis
                    plt.subplot(222)

                    plt.plot(timestamps[dp_c0b]-timestamps[0],mb[dp_c0b],'r-',drawstyle='steps-mid')
                    plt.plot(timestamps[dp_c0b]-timestamps[0],visb_ptr[dp_c0b,ch_plot],'g--',drawstyle='steps-mid')
                    plt.plot(timestamps[dp_c1b]-timestamps[0],mb[dp_c1b],'r-',drawstyle='steps-mid')
                    plt.plot(timestamps[dp_c1b]-timestamps[0],visb_ptr[dp_c1b,ch_plot],'g--',drawstyle='steps-mid')
                    #if fname in ['1551055211','1551037708','1579725085', '1580260015','1630519596']:
                    plt.plot(timestamps[dp_c2b]-timestamps[0],mb[dp_c2b],'r-',drawstyle='steps-mid')
                    plt.plot(timestamps[dp_c3b]-timestamps[0],mb[dp_c3b],'r-',drawstyle='steps-mid')
                    plt.plot(timestamps[dp_c4b]-timestamps[0],mb[dp_c4b],'r-',drawstyle='steps-mid')
                    plt.plot(timestamps[dp_c2b]-timestamps[0],visb_ptr[dp_c2b,ch_plot],'g--',drawstyle='steps-mid')
                    plt.plot(timestamps[dp_c3b]-timestamps[0],visb_ptr[dp_c3b,ch_plot],'g--',drawstyle='steps-mid')
                    plt.plot(timestamps[dp_c4b]-timestamps[0],visb_ptr[dp_c4b,ch_plot],'g--',drawstyle='steps-mid')
                    plt.ylim(np.nanmin(visb_ptr[dp_c0b,ch_plot])-20,np.nanmax(visb_ptr[dp_c1b,ch_plot])+30)
                    #plt.xlabel('time (s)')
                    #plt.ylabel('raw signal')
                    plt.legend(['model','raw signal'],ncol=2)
                    plt.title('Fitting curve of track-II at '+ str(int(freqs[ch_plot]/1e6)) +' MHz')
                    plt.subplot(224)

                    plt.plot([timestamps[dp_c0b][0]-timestamps[0],timestamps[dp_c0b][-1]-timestamps[0]],[0,0],'k-.')
                    plt.plot(timestamps[dp_c0b]-timestamps[0],resib[dp_c0b],'.-')
                    plt.plot(timestamps[dp_c1b]-timestamps[0],resib[dp_c1b],'.-')
                    #if fname in ['1551055211','1551037708','1579725085', '1580260015','1630519596']:
                    plt.plot(timestamps[dp_c2b]-timestamps[0],resib[dp_c2b],'.-')
                    plt.plot(timestamps[dp_c3b]-timestamps[0],resib[dp_c3b],'.-')
                    plt.plot(timestamps[dp_c4b]-timestamps[0],resib[dp_c4b],'.-')
                    plt.plot([timestamps[dp_c0b][0]-timestamps[0],timestamps[dp_c4b][-1]-timestamps[0]],[0,0],'k-.')
                    plt.xlabel('time (s)')
                    #plt.ylabel('residual (K)')
                    #if fname in ['1551055211','1551037708','1579725085', '1580260015','1630519596']:
                    plt.legend(['zero line','track c0', 'track c1', 'track c2', 'track c3','track c4'],fontsize=12, ncol=3)
                    #plt.title('calibrator Part II of '+str(fname)+', '+str(recv)+ ' @'+ str(round(freqs[ch_plot]/1e9,3)) +' GHz')
                    plt.ylim(-0.19,0.19)
                    plt.title('$T_{res}$ of the fitted curve above')
                    plt.savefig(output_file+'F_caliB_'+fname+'_'+recv+'_ch'+str(ch_plot)+'.png', bbox_inches='tight')
                    #plt.show()

                print (Tnd_ref)
                print (Tnda,Tndb,(Tnda+Tndb)/2.)
               
                ####data need to storage######
                Tnd_ref_list[ch_plot]=Tnd_ref
                Tel_map[:,ch_plot]=Tel
                
                if Tnda>0:
                    Tnda_list[ch_plot]=Tnda
                    NRMSE1_list[ch_plot]=NRMSE1
                if Tndb>0:
                    Tndb_list[ch_plot]=Tndb
                    NRMSE2_list[ch_plot]=NRMSE2
                    
                gta_param[ch_plot]=gta
                gtb_param[ch_plot]=gtb

                gain_map[dp_ca,ch_plot]=ga[dp_ca]
                gain_map[dp_cb,ch_plot]=gb[dp_cb]

                calT_tra=visa_ptr[:,ch_plot]/ga
                calT_trb=visb_ptr[:,ch_plot]/gb
                if Tnda>0:
                    assert((abs(calT_tra[dp_ca]-ma[dp_ca]/ga[dp_ca]-resia[dp_ca])<1e-10).all()==True)
                if Tndb>0:
                    assert((abs(calT_trb[dp_cb]-mb[dp_cb]/gb[dp_cb]-resib[dp_cb])<1e-10).all()==True)
                T_map[dp_ca,ch_plot]=calT_tra[dp_ca]
                T_map[dp_cb,ch_plot]=calT_trb[dp_cb]

                Tresi_map[dp_ca,ch_plot]=resia[dp_ca]
                Tresi_map[dp_cb,ch_plot]=resib[dp_cb]

                sma_param[ch_plot]=sma
                smb_param[ch_plot]=smb

                print (ch_plot, Tnd_ref, Tnda, Tndb)
                
                if target=='PKS1934-638':
                    ## extra flagger for weak calibrator

                    Tnda_raw,Tndb_raw=Tnda,Tndb

                    v_min1,v_max1=np.ma.min(visa_ptr[nd_0,ch_plot]),np.ma.max(visa_ptr[nd_0,ch_plot])
                    dp_list1=[dp_c1a,dp_c0a,dp_c2a,dp_c3a,dp_c4a]
                    vis_plot1=visa_ptr[:,ch_plot].copy()
                    vis_plot1.mask[nd_1x]=True

                    v_min2,v_max2=np.ma.min(visb_ptr[nd_0,ch_plot]),np.ma.max(visb_ptr[nd_0,ch_plot])
                    dp_list2=[dp_c1b,dp_c0b,dp_c2b,dp_c3b,dp_c4b]
                    vis_plot2=visb_ptr[:,ch_plot].copy()
                    vis_plot2.mask[nd_1x]=True

                    dp_flist1=[dp_c0a,dp_c1a1,dp_c2a,dp_c1a2,dp_c3a,dp_c1a3,dp_c4a]
                    mdn1_list,std1_list=kr.cal_std_list(dp_flist1,vis_plot1)
                    print (np.ma.mean(std1_list),std1_list)

                    dp_flist2=[dp_c0b,dp_c1b1,dp_c2b,dp_c1b2,dp_c3b,dp_c1b3,dp_c4b]
                    mdn2_list,std2_list=kr.cal_std_list(dp_flist2,vis_plot2)
                    print (np.ma.mean(std2_list),std2_list)

                    m_lim1=kr.set_vis_lim(dp_flist1,vis_plot1)
                    m_lim2=kr.set_vis_lim(dp_flist2,vis_plot2)

                    visa_ptr2=kr.filter_median(dp_flist1,mdn1_list,m_lim1,visa_ptr,ch_plot)
                    visb_ptr2=kr.filter_median(dp_flist2,mdn2_list,m_lim2,visb_ptr,ch_plot)

                    std_list=list(std1_list)+list(std2_list)
                    std_lim=np.min(std_list)
                    print (std_lim)

                    sigma_cut=2.5
                    visa_ptr2=kr.filter_std(dp_flist1,std1_list,visa_ptr2,std_lim,sigma_cut,ch_plot)
                    visb_ptr2=kr.filter_std(dp_flist2,std2_list,visb_ptr2,std_lim,sigma_cut,ch_plot)

                    #v_min1,v_max1=np.ma.min(visa_ptr2[nd_0,ch_plot]),np.ma.max(visa_ptr2[nd_0,ch_plot])
                    dp_list1=[dp_c1a,dp_c0a,dp_c2a,dp_c3a,dp_c4a]
                    vis_plot1=visa_ptr2[:,ch_plot].copy()
                    vis_plot1.mask[nd_1x]=True

                    #v_min2,v_max2=np.ma.min(visb_ptr2[nd_0,ch_plot]),np.ma.max(visb_ptr2[nd_0,ch_plot])
                    dp_list2=[dp_c1b,dp_c0b,dp_c2b,dp_c3b,dp_c4b]
                    vis_plot2=visb_ptr2[:,ch_plot].copy()
                    vis_plot2.mask[nd_1x]=True

                    track1_on_count,track1_off_count,track2_on_count,track2_off_count=3,4,3,4
                    for i in [dp_c1a1,dp_c1a2,dp_c1a3]:
                        if (vis_plot1[i].mask==True).all():
                            track1_on_count-=1
                    for i in [dp_c0a,dp_c2a,dp_c3a,dp_c4a]:
                        if (vis_plot1[i].mask==True).all():
                            track1_off_count-=1
                    for i in [dp_c1b1,dp_c1b2,dp_c1b3]:
                        if (vis_plot2[i].mask==True).all():
                            track2_on_count-=1
                    for i in [dp_c0b,dp_c2b,dp_c3b,dp_c4b]:
                        if (vis_plot2[i].mask==True).all():
                            track2_off_count-=1

                    print(track1_on_count,track1_off_count,track2_on_count,track2_off_count)

                    if track1_on_count>1 and track1_off_count>1:
                        Trec0=Trec_list[ch_plot]
                        print (Trec0)
                        eta_p0=1.0
                        func_sm_param0=[Trec0]
                        func_gt_param0=[ga0,0,0,0,0]
                        print (Tnd_ref)
                        ####fitting
                        instru_pa=ks.solve_params0_v3(timestamps, visa_ptr2, ch_plot, nd_ratio, T_ptr, eta_p0, Tnd_ref, Tnd_std, Tel, Tgal, func_gt_param0, func_sm_param0, nd_0, nd_1x)

                        ####get fitting result
                        Tnda=instru_pa[0]
                        eta_pa=instru_pa[1]
                        sma=instru_pa[2]
                        gta=instru_pa[3:]

                        print (Tnda, eta_pa, sma, gta)
                    else:
                        Tnda, eta_pa, sma, gta=0,0,0,0
                        print ('# not suitable for fitting')

                    if track2_on_count>1 and track2_off_count>1:

                        print (Trec0)
                        eta_p0=1.0
                        func_sm_param0=[Trec0]
                        func_gt_param0=[gb0,0,0,0,0]

                        print (Tnd_ref)

                        ####fitting######
                        instru_pb=ks.solve_params0_v3(timestamps, visb_ptr2, ch_plot, nd_ratio, T_ptr, eta_p0, Tnd_ref, Tnd_std, Tel, Tgal,
                                              func_gt_param0, func_sm_param0, nd_0, nd_1x)

                        ######get fitting result#####
                        Tndb=instru_pb[0]
                        eta_pb=instru_pb[1]
                        smb=instru_pb[2]
                        gtb=instru_pb[3:]

                        print (Tndb, eta_pb, smb, gtb)
                    else:
                        Tndb, eta_pb, smb, gtb=0,0,0,0
                        print ('# not suitable for fitting')


                    m=ks.calc_total_model_v3(timestamps, nd_ratio, T_ptr, eta_pa, Tnda, Tel, Tgal, gta, sma, nd_0, nd_1x)
                    ma=np.ma.array(m,mask=visa_ptr[:,ch_plot].mask)
                    m=ks.calc_total_model_v3(timestamps, nd_ratio, T_ptr, eta_pb, Tndb, Tel, Tgal, gtb, smb, nd_0, nd_1x)
                    mb=np.ma.array(m,mask=visb_ptr[:,ch_plot].mask)
                    ga=ks.func_gt(timestamps,gta)
                    resia=(visa_ptr2[:,ch_plot]-ma)/ga
                    gb=ks.func_gt(timestamps,gtb)
                    resib=(visb_ptr2[:,ch_plot]-mb)/gb
                    
                    NRMSE11=ks.cal_NRMSE(ma/ga,resia)
                    NRMSE22=ks.cal_NRMSE(mb/gb,resib)

                    if ch_plot==800:
                        ##show model and raw vis
                        plt.figure(figsize=(20,8))
                        plt.subplots_adjust(wspace=0.1,hspace=0.25)
                        plt.subplot(221)

                        plt.plot(timestamps[dp_c0a]-timestamps[0],ma[dp_c0a],'r-',drawstyle='steps-mid')
                        plt.plot(timestamps[dp_c0a]-timestamps[0],visa_ptr2[dp_c0a,ch_plot],'g--',drawstyle='steps-mid')
                        plt.plot(timestamps[dp_c1a]-timestamps[0],ma[dp_c1a],'r-',drawstyle='steps-mid')
                        plt.plot(timestamps[dp_c1a]-timestamps[0],visa_ptr2[dp_c1a,ch_plot],'g--',drawstyle='steps-mid')
                        #if fname in ['1551055211','1551037708','1630519596']:
                        plt.plot(timestamps[dp_c2a]-timestamps[0],ma[dp_c2a],'r-',drawstyle='steps-mid')
                        plt.plot(timestamps[dp_c3a]-timestamps[0],ma[dp_c3a],'r-',drawstyle='steps-mid')
                        plt.plot(timestamps[dp_c4a]-timestamps[0],ma[dp_c4a],'r-',drawstyle='steps-mid')
                        plt.plot(timestamps[dp_c2a]-timestamps[0],visa_ptr2[dp_c2a,ch_plot],'g--',drawstyle='steps-mid')
                        plt.plot(timestamps[dp_c3a]-timestamps[0],visa_ptr2[dp_c3a,ch_plot],'g--',drawstyle='steps-mid')
                        plt.plot(timestamps[dp_c4a]-timestamps[0],visa_ptr2[dp_c4a,ch_plot],'g--',drawstyle='steps-mid')
                        #plt.xlabel('time (s)')
                        plt.ylabel('raw signal')
                        plt.legend(['model','raw signal'],ncol=2)
                        plt.ylim(np.nanmin(visa_ptr[dp_c0a,ch_plot])-20,np.nanmax(visa_ptr[dp_c1a,ch_plot])+30)
                        plt.title('Fitting curve2 of track-I at '+ str(int(freqs[ch_plot]/1e6)) +' MHz, '+str(fname)+', '+str(recv))
                        if fname=='1551055211':
                            plt.title('Fitting curve of track-I at '+ str(int(freqs[ch_plot]/1e6)) +' MHz, '+str(recv)+' of obs190225')
                        plt.subplot(223)

                        plt.plot([timestamps[dp_c0a][0]-timestamps[0],timestamps[dp_c0a][-1]-timestamps[0]],[0,0],'k-.')
                        plt.plot(timestamps[dp_c0a]-timestamps[0],resia[dp_c0a],'.-')
                        plt.plot(timestamps[dp_c1a]-timestamps[0],resia[dp_c1a],'.-')
                        #if fname in ['1551055211','1551037708','1630519596']:
                        plt.plot(timestamps[dp_c2a]-timestamps[0],resia[dp_c2a],'.-')
                        plt.plot(timestamps[dp_c3a]-timestamps[0],resia[dp_c3a],'.-')
                        plt.plot(timestamps[dp_c4a]-timestamps[0],resia[dp_c4a],'.-')
                        plt.plot([timestamps[dp_c0a][0]-timestamps[0],timestamps[dp_c4a][-1]-timestamps[0]],[0,0],'k-.')
                        plt.xlabel('time (s)')
                        plt.ylabel('Tres (K)')
                        plt.ylim(-0.15,0.15)
                        #plt.legend(['line of T=0 K','track calibrator outskirt', 'track calibrator center'],ncol=3)
                        #if fname in ['1551055211','1551037708','1579725085', '1580260015','1630519596']:
                        plt.legend(['zero line','track c0', 'track c1', 'track c2', 'track c3','track c4'], fontsize=12, ncol=3)
                        plt.ylim(-0.19,0.19)
                        plt.title('$T_{res}$ of the fitted curve above')
                        #plt.savefig('caliA_ch'+str(ch_plot)+'.png')
                        ##show model and raw vis
                        plt.subplot(222)

                        plt.plot(timestamps[dp_c0b]-timestamps[0],mb[dp_c0b],'r-',drawstyle='steps-mid')
                        plt.plot(timestamps[dp_c0b]-timestamps[0],visb_ptr2[dp_c0b,ch_plot],'g--',drawstyle='steps-mid')
                        plt.plot(timestamps[dp_c1b]-timestamps[0],mb[dp_c1b],'r-',drawstyle='steps-mid')
                        plt.plot(timestamps[dp_c1b]-timestamps[0],visb_ptr2[dp_c1b,ch_plot],'g--',drawstyle='steps-mid')
                        #if fname in ['1551055211','1551037708','1579725085', '1580260015','1630519596']:
                        plt.plot(timestamps[dp_c2b]-timestamps[0],mb[dp_c2b],'r-',drawstyle='steps-mid')
                        plt.plot(timestamps[dp_c3b]-timestamps[0],mb[dp_c3b],'r-',drawstyle='steps-mid')
                        plt.plot(timestamps[dp_c4b]-timestamps[0],mb[dp_c4b],'r-',drawstyle='steps-mid')
                        plt.plot(timestamps[dp_c2b]-timestamps[0],visb_ptr2[dp_c2b,ch_plot],'g--',drawstyle='steps-mid')
                        plt.plot(timestamps[dp_c3b]-timestamps[0],visb_ptr2[dp_c3b,ch_plot],'g--',drawstyle='steps-mid')
                        plt.plot(timestamps[dp_c4b]-timestamps[0],visb_ptr2[dp_c4b,ch_plot],'g--',drawstyle='steps-mid')
                        plt.ylim(np.nanmin(visb_ptr2[dp_c0b,ch_plot])-20,np.nanmax(visb_ptr2[dp_c1b,ch_plot])+30)
                        #plt.xlabel('time (s)')
                        #plt.ylabel('raw signal')
                        plt.legend(['model','raw signal'],ncol=2)
                        plt.title('Fitting curve of track-II at '+ str(int(freqs[ch_plot]/1e6)) +' MHz')
                        plt.subplot(224)

                        plt.plot([timestamps[dp_c0b][0]-timestamps[0],timestamps[dp_c0b][-1]-timestamps[0]],[0,0],'k-.')
                        plt.plot(timestamps[dp_c0b]-timestamps[0],resib[dp_c0b],'.-')
                        plt.plot(timestamps[dp_c1b]-timestamps[0],resib[dp_c1b],'.-')
                        #if fname in ['1551055211','1551037708','1579725085', '1580260015','1630519596']:
                        plt.plot(timestamps[dp_c2b]-timestamps[0],resib[dp_c2b],'.-')
                        plt.plot(timestamps[dp_c3b]-timestamps[0],resib[dp_c3b],'.-')
                        plt.plot(timestamps[dp_c4b]-timestamps[0],resib[dp_c4b],'.-')
                        plt.plot([timestamps[dp_c0b][0]-timestamps[0],timestamps[dp_c4b][-1]-timestamps[0]],[0,0],'k-.')
                        plt.xlabel('time (s)')
                        #plt.ylabel('residual (K)')
                        #if fname in ['1551055211','1551037708','1579725085', '1580260015','1630519596']:
                        plt.legend(['zero line','track c0', 'track c1', 'track c2', 'track c3','track c4'],fontsize=12, ncol=3)
                        #plt.title('calibrator Part II of '+str(fname)+', '+str(recv)+ ' @'+ str(round(freqs[ch_plot]/1e9,3)) +' GHz')
                        plt.ylim(-0.19,0.19)
                        plt.title('$T_{res}$ of the fitted curve above')
                        plt.savefig(output_file+'F_caliB_'+fname+'_'+recv+'_ch'+str(ch_plot)+'_2.png', bbox_inches='tight')
                        #plt.show()

                    print(track1_on_count,track1_off_count,track2_on_count,track2_off_count)

                    print ('Tnd_ref: '+str(Tnd_ref))
                    print ('pre-fit: '+str(Tnda_raw)+' '+str(Tndb_raw))
                    if Tnda>0:
                        print ('Tnda: '+str(Tnda))
                    if Tndb>0:
                        print ('Tndb: '+str(Tndb))
                    if Tnda>0 and Tndb>0:
                        print ('Tnd mean: '+str((Tnda+Tndb)/2.))

                    ####data need to storage######
                    
                    if Tnda>0:
                        Tnda_list2[ch_plot]=Tnda
                        NRMSE11_list[ch_plot]=NRMSE11
                    if Tndb>0:
                        Tndb_list2[ch_plot]=Tndb
                        NRMSE22_list[ch_plot]=NRMSE22
                    
                    gta_param2[ch_plot]=gta
                    gtb_param2[ch_plot]=gtb

                    gain_map2[dp_ca,ch_plot]=ga[dp_ca]
                    gain_map2[dp_cb,ch_plot]=gb[dp_cb]

                    calT_tra=visa_ptr2[:,ch_plot]/ga
                    calT_trb=visb_ptr2[:,ch_plot]/gb
                    if Tnda>0:
                        assert((abs(calT_tra[dp_ca]-ma[dp_ca]/ga[dp_ca]-resia[dp_ca])<1e-10).all()==True)
                    if Tndb>0:
                        assert((abs(calT_trb[dp_cb]-mb[dp_cb]/gb[dp_cb]-resib[dp_cb])<1e-10).all()==True)
                    T_map2[dp_ca,ch_plot]=calT_tra[dp_ca]
                    T_map2[dp_cb,ch_plot]=calT_trb[dp_cb]

                    Tresi_map2[dp_ca,ch_plot]=resia[dp_ca]
                    Tresi_map2[dp_cb,ch_plot]=resib[dp_cb]

                    sma_param2[ch_plot]=sma
                    smb_param2[ch_plot]=smb

                    print (ch_plot, Tnd_ref, Tnda, Tndb)

                if ch_plot%50==0:
                    ####save data####
                    d['timestamps']=timestamps
                    d['nd_0']=nd_0
                    d['ra']=ra
                    d['dec']=dec
                    d['Tnd_ref_list']=Tnd_ref_list
                    d['Tel_map']=Tel_map

                    d['T_map']=T_map
                    d['Tresi_map']=Tresi_map
                    d['gain_map']=gain_map
                    d['Tnda_list']=Tnda_list
                    d['Tndb_list']=Tndb_list
                    d['gta_param']=gta_param
                    d['gtb_param']=gtb_param
                    d['sma_param']=sma_param
                    d['smb_param']=smb_param
                    d['NRMSE1_list']=NRMSE1_list
                    d['NRMSE2_list']=NRMSE2_list

                    if target=='PKS1934-638':
                        d['T_map2']=T_map2
                        d['Tresi_map2']=Tresi_map2
                        d['gain_map2']=gain_map2
                        d['Tnda_list2']=Tnda_list2
                        d['Tndb_list2']=Tndb_list2
                        d['gta_param2']=gta_param2
                        d['gtb_param2']=gtb_param2
                        d['sma_param2']=sma_param2
                        d['smb_param2']=smb_param2
                        d['NRMSE11_list']=NRMSE11_list
                        d['NRMSE22_list']=NRMSE22_list
                    fs=open(output_file+str(fname)+'_'+str(recv)+'_level2_data','wb')
                    pickle.dump(d,fs,protocol=2)
                    fs.close()

                    d2['Tnd_ref_list']=Tnd_ref_list
                    d2['Tnda_list']=Tnda_list
                    d2['Tndb_list']=Tndb_list
                    d2['NRMSE1_list']=NRMSE1_list
                    d2['NRMSE2_list']=NRMSE2_list
                    if target=='PKS1934-638':
                        d2['Tnda_list2']=Tnda_list2
                        d2['Tndb_list2']=Tndb_list2
                        d2['NRMSE11_list']=NRMSE11_list
                        d2['NRMSE22_list']=NRMSE22_list
                    fs=open(output_file+str(fname)+'_'+str(recv)+'_level2_Tnd_data','wb')
                    pickle.dump(d2,fs,protocol=2)
                    fs.close()
                    print ('tmp data saved @ch'+str(ch_plot))
                print ('***channel'+ str(ch_plot) +' finished')    
                count_ch+=1
                
            except(Exception):
                print ('***channel'+ str(ch_plot) +' failed')
        else:
            print ('***channel'+ str(ch_plot) +' skipped')        


    ####save data####
    d['timestamps']=timestamps
    d['nd_0']=nd_0
    d['ra']=ra
    d['dec']=dec
    d['Tnd_ref_list']=Tnd_ref_list
    d['Tel_map']=Tel_map
    
    d['T_map']=T_map
    d['Tresi_map']=Tresi_map
    d['gain_map']=gain_map
    d['Tnda_list']=Tnda_list
    d['Tndb_list']=Tndb_list
    d['gta_param']=gta_param
    d['gtb_param']=gtb_param
    d['sma_param']=sma_param
    d['smb_param']=smb_param
    d['NRMSE1_list']=NRMSE1_list
    d['NRMSE2_list']=NRMSE2_list
    
    if target=='PKS1934-638':
        d['T_map2']=T_map2
        d['Tresi_map2']=Tresi_map2
        d['gain_map2']=gain_map2
        d['Tnda_list2']=Tnda_list2
        d['Tndb_list2']=Tndb_list2
        d['gta_param2']=gta_param2
        d['gtb_param2']=gtb_param2
        d['sma_param2']=sma_param2
        d['smb_param2']=smb_param2
        d['NRMSE11_list']=NRMSE11_list
        d['NRMSE22_list']=NRMSE22_list
    fs=open(output_file+str(fname)+'_'+str(recv)+'_level2_data','wb')
    pickle.dump(d,fs,protocol=2)
    fs.close()

    
    
    ###begin: reference plot only###
    NRMSE_lim=0.004

    Tnda_list=kio.mask_None(Tnda_list)
    Tndb_list=kio.mask_None(Tndb_list)
    NRMSE1_list=kio.mask_None(NRMSE1_list)
    NRMSE2_list=kio.mask_None(NRMSE2_list)
    if target=='PKS1934-638':
        Tnda_list2=kio.mask_None(Tnda_list2)
        Tndb_list2=kio.mask_None(Tndb_list2)
        NRMSE11_list=kio.mask_None(NRMSE11_list)
        NRMSE22_list=kio.mask_None(NRMSE22_list)

    plt.figure(figsize=(16,5))
    plt.plot(Tnda_list,'.')
    plt.plot(Tndb_list,'.')
    plt.plot(Tnd_ref_list, color='grey',zorder=0)
    if target=='PKS1934-638':
        plt.plot(Tnda_list2,'x',zorder=1)
        plt.plot(Tndb_list2,'x',zorder=2)

    plt.legend(['Tnda','Tndb','Tnd_ref','Tnda2','Tndb2'],ncol=2)
    plt.xlabel('channel')
    plt.ylabel('Tnd (K)')
    plt.title('Tnd result '+fname+', '+recv)
    plt.xlim(550,1050)
    plt.savefig(output_file+'Tnd_all_'+str(fname)+'_'+str(recv)+'_l.png', bbox_inches='tight')
    
    plt.figure(figsize=(16,5))
    plt.plot(Tnda_list,'.')
    plt.plot(Tndb_list,'.')
    plt.plot(Tnd_ref_list, color='grey',zorder=0)
    if target=='PKS1934-638':
        plt.plot(Tnda_list2,'x',zorder=1)
        plt.plot(Tndb_list2,'x',zorder=2)

    plt.legend(['Tnda','Tndb','Tnd_ref','Tnda2','Tndb2'],ncol=2)
    plt.xlabel('channel')
    plt.ylabel('Tnd (K)')
    plt.title('Tnd result '+fname+', '+recv)
    plt.xlim(2150,3100)
    plt.savefig(output_file+'Tnd_all_'+str(fname)+'_'+str(recv)+'_h.png', bbox_inches='tight')
    
    #plt.show()

    sigma_local=3.

    d_list=kf.curve_filter_ma_array(range(4096),Tnda_list,sigma=sigma_local)
    Tnda_list.mask[d_list]=True
    d_list=kf.curve_filter_ma_array(range(4096),Tndb_list,sigma=sigma_local)
    Tndb_list.mask[d_list]=True
    if target=='PKS1934-638':
        d_list=kf.curve_filter_ma_array(range(4096),Tnda_list2,sigma=sigma_local)
        Tnda_list2.mask[d_list]=True
        d_list=kf.curve_filter_ma_array(range(4096),Tndb_list2,sigma=sigma_local)
        Tndb_list2.mask[d_list]=True

    Tnda_list.mask[NRMSE1_list>NRMSE_lim]=True
    Tndb_list.mask[NRMSE2_list>NRMSE_lim]=True
    
    Tnd_lim=10
    Tnda_list.mask[Tnda_list<Tnd_lim]=True
    Tndb_list.mask[Tndb_list<Tnd_lim]=True
    if target=='PKS1934-638':
        Tnda_list2.mask[NRMSE11_list>NRMSE_lim]=True
        Tndb_list2.mask[NRMSE22_list>NRMSE_lim]=True
        Tnda_list2.mask[Tnda_list2<Tnd_lim]=True
        Tndb_list2.mask[Tndb_list2<Tnd_lim]=True
        
    #if target=='PKS1934-638':
        Tnda_final=np.ma.array(np.zeros(4096),mask=True)
        Tndb_final=np.ma.array(np.zeros(4096),mask=True)

        mask1=Tnda_list.mask
        mask11=Tnda_list2.mask
        Tnda_final[(~mask1) & (~mask11) & (NRMSE1_list<= NRMSE11_list)]=Tnda_list[(~mask1) & (~mask11) & (NRMSE1_list<= NRMSE11_list)]
        Tnda_final[(~mask1) & (~mask11) & (NRMSE1_list> NRMSE11_list)]=Tnda_list2[(~mask1) & (~mask11) & (NRMSE1_list> NRMSE11_list)]
        Tnda_final[(~mask1) & (mask11)]=Tnda_list[(~mask1) & (mask11)]
        Tnda_final[(mask1) & (~mask11)]=Tnda_list2[(mask1) & (~mask11)]

        mask2=Tndb_list.mask
        mask22=Tndb_list2.mask
        Tndb_final[(~mask2) & (~mask22) & (NRMSE2_list<= NRMSE22_list)]=Tndb_list[(~mask2) & (~mask22) & (NRMSE2_list<= NRMSE22_list)]
        Tndb_final[(~mask2) & (~mask22) & (NRMSE2_list> NRMSE22_list)]=Tndb_list2[(~mask2) & (~mask22) & (NRMSE2_list> NRMSE22_list)]
        Tndb_final[(~mask2) & (mask22)]=Tndb_list[(~mask2) & (mask22)]
        Tndb_final[(mask2) & (~mask22)]=Tndb_list2[(mask2) & (~mask22)]
    else:
        Tnda_final=Tnda_list.copy()
        Tndb_final=Tndb_list.copy()

    Tnd_final=np.ma.mean([Tnda_final,Tndb_final],axis=0)
    #Tnd_final.mask[abs(Tnda_final-Tndb_final)>1]=True #not what you think, below is the right version
    Tnd_final.mask[np.ma.where(np.ma.abs(Tnda_final-Tndb_final)>1)[0]]=True
    
    
    plt.figure(figsize=(16,5))
    plt.plot(Tnda_final,'.')
    plt.plot(Tndb_final,'.')
    plt.plot(Tnd_final,'x')
    plt.plot(Tnd_ref_list, color='grey',zorder=0)
    plt.legend(['Tnda_final','Tndb_final','Tnd_final','Tnd_ref'],ncol=2)
    plt.xlabel('channel')
    plt.ylabel('Tnd (K)')
    plt.title('Tnd result '+fname+', '+recv)
    plt.xlim(550,1050)
    plt.savefig(output_file+'Tnd_final_pre_'+str(fname)+'_'+str(recv)+'_l.png', bbox_inches='tight')
    
    plt.figure(figsize=(16,5))
    plt.plot(Tnda_final,'.')
    plt.plot(Tndb_final,'.')
    plt.plot(Tnd_final,'x')
    plt.plot(Tnd_ref_list, color='grey',zorder=0)
    plt.legend(['Tnda_final','Tndb_final','Tnd_final','Tnd_ref'],ncol=2)
    plt.xlabel('channel')
    plt.ylabel('Tnd (K)')
    plt.title('Tnd result '+fname+', '+recv)
    plt.xlim(2150,3100)
    plt.savefig(output_file+'Tnd_final_pre_'+str(fname)+'_'+str(recv)+'_h.png', bbox_inches='tight')
    
    #plt.show()
    
    sigma_local2=2.5

    d_list=kf.curve_filter_ma_array(range(4096),Tnd_final,sigma=sigma_local2)
    Tnd_final.mask[d_list]=True
    
    d2['Tnd_final']=Tnd_final #new added
    d2['Tnd_ref_list']=Tnd_ref_list
    d2['Tnda_list']=Tnda_list
    d2['Tndb_list']=Tndb_list
    d2['NRMSE1_list']=NRMSE1_list
    d2['NRMSE2_list']=NRMSE2_list
    if target=='PKS1934-638':
        d2['Tnda_list2']=Tnda_list2
        d2['Tndb_list2']=Tndb_list2
        d2['NRMSE11_list']=NRMSE11_list
        d2['NRMSE22_list']=NRMSE22_list
    fs=open(output_file+str(fname)+'_'+str(recv)+'_level2_Tnd_data','wb')
    pickle.dump(d2,fs,protocol=2)
    fs.close()
    
   
    plt.figure(figsize=(16,5))
    plt.plot(Tnd_final,'x')
    plt.plot(Tnd_ref_list, color='grey',zorder=0)
    plt.legend(['Tnd_final','Tnd_ref'],ncol=2)
    plt.xlabel('channel')
    plt.ylabel('Tnd (K)')
    plt.title('Tnd result '+fname+', '+recv)
    plt.xlim(550,1050)
    plt.savefig(output_file+'Tnd_final_'+str(fname)+'_'+str(recv)+'_l.png', bbox_inches='tight')
    plt.show()
    
    plt.figure(figsize=(16,5))
    plt.plot(Tnd_final,'x')
    plt.plot(Tnd_ref_list, color='grey',zorder=0)
    plt.legend(['Tnd_final','Tnd_ref'],ncol=2)
    plt.xlabel('channel')
    plt.ylabel('Tnd (K)')
    plt.title('Tnd result '+fname+', '+recv)
    plt.xlim(2150,3100)
    plt.savefig(output_file+'Tnd_final_'+str(fname)+'_'+str(recv)+'_h.png', bbox_inches='tight')
    plt.show()
    
    ###end: referncer plot only###
    print ('total calibrated channel count:'+str(count_ch))
    print ('***Level2 calibartion for '+fname+' '+ant+pol)    
    print ('end @ ' + time.asctime(time.localtime(time.time())) +'#')
else:
    print ('***BAD ANT, Level2 calibartion for '+fname+' '+ant+pol+' skipped')
    print ('end @ ' + time.asctime(time.localtime(time.time())) +'#')
print ('Hakuna Matata') 
################################################
#### Jingying Wang <astro.jywang@gmail.com> ####
################################################
