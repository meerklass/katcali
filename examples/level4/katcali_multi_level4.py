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

p_radec=np.loadtxt('radio_source.txt')

output_file='./level4_output/'

def cal_map_I(map_h, map_v):    
    assert(np.shape(map_h)==np.shape(map_v))
    map=(map_h+map_v)/2.   
    print 'min value h,v:',
    print np.ma.min(map_h),np.ma.min(map_v)
    assert(np.ma.min(map_h)>0)
    assert(np.ma.min(map_v)>0)
    diff_ratio=np.ma.exp(abs(np.ma.log(map_h)-np.ma.log(map_v))) #difference betwwen map h,v; 
                                                        #function to make sure maps are exchangable
    return map, diff_ratio



fname,ant=sys.argv[1],sys.argv[2]

print fname,ant

try:
    dh=pickle.load(open('/idia/projects/hi_im/raw_vis/katcali_output/level3_output/'+str(fname)+'_'+str(ant)+'h_level3_data'))
    dv=pickle.load(open('/idia/projects/hi_im/raw_vis/katcali_output/level3_output/'+str(fname)+'_'+str(ant)+'v_level3_data'))

    assert((dh['ra']==dv['ra']).all()==True)
    assert((dh['dec']==dv['dec']).all()==True)
    assert((dh['timestamps']==dv['timestamps']).all()==True)
    assert((dh['Tgal_map']==dv['Tgal_map']).all()==True)
    ra=dh['ra']
    dec=dh['dec']
    timestamps=dh['timestamps']
    Tgal_map=dh['Tgal_map']

    gain_map_h=dh['gain_map']
    T_map_h=dh['T_map']
    Tresi_map_h=dh['Tresi_map']
    Tsm_map_h=dh['Tsm_map']
    mask_nd_h=dh['mask_nd_s0'] #True is diode ON
    Tel_map_h=dh['Tel_map']
    Tsky_map_h=T_map_h-Tel_map_h-Tsm_map_h

    gain_map_v=dv['gain_map']
    T_map_v=dv['T_map']
    Tresi_map_v=dv['Tresi_map']
    Tsm_map_v=dv['Tsm_map']
    mask_nd_v=dv['mask_nd_s0'] #True is diode ON
    Tel_map_v=dv['Tel_map']
    Tsky_map_v=T_map_v-Tel_map_v-Tsm_map_v


    ch_e=3101


    assert((T_map_h.mask==Tresi_map_h.mask).all()==True)
    assert((T_map_v.mask==Tresi_map_v.mask).all()==True)
    div=T_map_h/T_map_v
    mask_clean=div.mask.copy()
    mask_clean_backup=mask_clean.copy()

    l=[T_map_h,T_map_v,Tresi_map_h,Tresi_map_v,Tsm_map_h,Tsm_map_v,gain_map_h,gain_map_v]
    l_str=['T_map_h','T_map_v','Tresi_map_h','Tresi_map_v','Tsm_map_h','Tsm_map_v','gain_map_h','gain_map_v']

    niter=6

    plt.figure(figsize=(10,len(l)*niter))
    plt.subplots_adjust (wspace=0.2, hspace=0.6) 
    for i in range(niter):
        print '> The iteration '+str(i+1)+' is in progress...'
        mask_clean2=mask_clean.copy()

        del_point=[]

        for l_i in range(len(l)):
            a=l[l_i]
            plt.subplot(len(l)/2*niter,2,i*len(l)+l_i+1)
            a_ch=np.ma.mean(np.ma.array(a,mask=mask_clean),axis=0)
            plt.plot(a_ch,'g.')

            ax1=kf.curve_filter_ma_array(range(550,1051),a_ch[550:1051],sigma=4.)#550-1050 filter
            plt.plot(ax1,a_ch[ax1],'r.') #deleted data
            mask_clean2[:,ax1]=True ###apply to the mask
            del_point+=list(ax1)

            ax2=kf.curve_filter_ma_array(range(2150,3101),a_ch[2150:3101],sigma=4.)#2150-3100 filter
            plt.plot(ax2,a_ch[ax2],'m.') #deleted data
            mask_clean2[:,ax2]=True ###apply to the mask
            del_point+=list(ax2)

            ax3=kf.curve_filter_ma_array(range(550,3101),a_ch[550:3101],sigma=4.,k=1)
            plt.plot(ax3,a_ch[ax3],'x',c='orange') #deleted data
            mask_clean2[:,ax3]=True ###apply to the mask
            del_point+=list(ax3)

            plt.title(l_str[l_i]+', iter '+str(i+1))

        print 'masks before and after filter are same:',
        print (mask_clean==mask_clean2).all()
        mask_clean=mask_clean2.copy()

        print str(len(set(del_point))) + ' channels masked by '+str(len(del_point))+' times'
        if len(del_point)==0:
            print '***filter done after iter time '+str(i+1) 
            break
    plt.savefig(output_file+fname+'_'+ant+'_ch_filter_iter'+str(i+1)+'.png')
    #plt.show()

    Tsky_map_h=np.ma.array(Tsky_map_h,mask=mask_clean)
    Tsky_map_v=np.ma.array(Tsky_map_v,mask=mask_clean)

    Tresi_map_h=np.ma.array(Tresi_map_h,mask=mask_clean)
    Tresi_map_v=np.ma.array(Tresi_map_v,mask=mask_clean)


    Tsky_map,Tsky_ratio=cal_map_I(Tsky_map_h, Tsky_map_v)
    Tresi_map=(Tresi_map_h+Tresi_map_v)/2.
    assert((Tsky_map.mask==Tresi_map.mask).all()==True)

    plt.figure(figsize=(12,4.5))
    plt.subplot(121)
    plt.imshow(Tsky_map,aspect='auto')
    plt.xlim(550,ch_e)
    plt.title(fname+'_'+ant+'_Tsky_map')
    plt.colorbar()
    plt.subplot(122)
    plt.imshow(Tresi_map,aspect='auto')
    plt.xlim(550,ch_e)
    plt.title(fname+'_'+ant+'_Tresi_map')
    plt.colorbar()
    plt.savefig(output_file+fname+'_'+ant+'_Inten.png')
    #plt.show()



    assert((Tsky_ratio.mask==mask_clean).all()==True)
    print np.ma.max(Tsky_ratio)

    plt.figure(figsize=(6,4.5))
    plt.imshow(Tsky_ratio,aspect='auto')
    plt.xlim(550,ch_e)
    plt.title(fname+'_'+ant+'_Tsky_ratio')
    plt.colorbar()
    plt.savefig(output_file+fname+'_'+ant+'_Tsky_ratio.png')
    #plt.show()

    bins=200

    plt.figure(figsize=(10,3))
    max=np.sum(Tsky_ratio.mask==False)/bins
    median=np.ma.median(Tsky_ratio)
    plt.hist(Tsky_ratio[~Tsky_ratio.mask],bins=bins,log=True)
    plt.plot([median,median],[0,max],'g-',lw=5)
    plt.xlabel('Tsky ratio')
    plt.ylabel('Count')
    plt.title(fname+'_'+ant+', bins'+str(bins))
    plt.savefig(output_file+fname+'_'+ant+'_Tsky_ratio_his.png')
    #plt.show()

    ch_plot_list=[600,700,800,900,1000,2200,2400,2600,2800,3000]
    plt.figure(figsize=(24,20))
    plt.subplots_adjust (wspace=0.0, hspace=0.2) 
    for i in range(len(ch_plot_list)):
        ch_plot1=ch_plot_list[i]
        p_data=Tsky_map[:,ch_plot1]
        plt.subplot(5,2,i+1)
        try:
            kv.plot_data(ra,dec, p_data,gsize=90)
        except(Exception):
            kv.plot_data(ra,dec, np.zeros_like(ra),gsize=90)
        if i>len(ch_plot_list)-3:
            plt.xlabel('R.A.')
        plt.ylabel('Dec')
        plt.title(fname+'_'+ant+'_ch'+str(ch_plot1))
    plt.savefig(output_file+fname+'_'+ant+'_Tsky_map.png')
    #plt.show()


    plt.figure(figsize=(24,20))
    plt.subplots_adjust (wspace=0.0, hspace=0.2) 
    for i in range(len(ch_plot_list)):
        ch_plot1=ch_plot_list[i]
        p_data=Tresi_map[:,ch_plot1]
        plt.subplot(5,2,i+1)
        try:
            kv.plot_data(ra,dec, p_data,gsize=90)
        except(Exception):
            kv.plot_data(ra,dec, np.zeros_like(ra),gsize=90)
        if i>len(ch_plot_list)-3:
            plt.xlabel('R.A.')
        plt.ylabel('Dec')
        plt.title(fname+'_'+ant+'_ch'+str(ch_plot1))
    plt.savefig(output_file+fname+'_'+ant+'_Tresi_map.png')
    #plt.show()


    plot_gsize=90
    plt.figure(figsize=(22,7))
    plt.subplots_adjust(wspace=0,hspace=.3)

    plt.subplot(221)
    p_data=np.ma.mean(Tsky_map,axis=1)
    plt.scatter(ra,dec, c=p_data, vmin=p_data.min(),vmax=p_data.max(), s=8)
    plt.ylabel('Dec (J2000)',fontsize=12)
    plt.title(fname+' '+ant+', fband Tsky_map')
    plt.colorbar()
    plt.subplot(222)
    kv.plot_data(ra,dec, p_data,gsize=plot_gsize)
    plt.plot(p_radec[:,0],p_radec[:,1],'mo')
    plt.ylabel('Dec (J2000)',fontsize=12)

    plt.subplot(223)
    p_data=np.ma.mean(Tresi_map,axis=1)
    plt.scatter(ra,dec, c=p_data, vmin=p_data.min(),vmax=p_data.max(), s=8)
    plt.ylabel('Dec (J2000)',fontsize=12)
    plt.title(fname+' '+ant+', fband Tresi_map')
    plt.colorbar()
    plt.subplot(224)
    kv.plot_data(ra,dec, p_data,gsize=plot_gsize)
    plt.plot(p_radec[:,0],p_radec[:,1],'mo')
    plt.ylabel('Dec (J2000)',fontsize=12)
    plt.savefig(output_file+fname+'_'+ant+'_fband_map.png')
    #plt.show()

    assert((Tsky_map.mask==mask_clean).all()==True)
    assert((Tresi_map.mask==mask_clean).all()==True)
    assert((Tsky_ratio.mask==mask_clean).all()==True)

    d={}
    d['Tsky_map']=Tsky_map
    d['Tresi_map']=Tresi_map
    d['Tsky_ratio']=Tsky_ratio

    d['timestamps']=timestamps
    d['ra']=ra
    d['dec']=dec
    fs=open(output_file+str(fname)+'_'+str(ant)+'_level4_data','wb')
    pickle.dump(d,fs,protocol=2)
    fs.close()

    d1={}

    d1['Inten_mask']=mask_clean
    fs=open(output_file+str(fname)+'_'+str(ant)+'_level4_mask','wb')
    pickle.dump(d1,fs,protocol=2)
    fs.close()

except IOError:
    print '#### no pre data exists.'


print 'end @ ' + time.asctime(time.localtime(time.time())) +'#'
print 'Hakuna Matata' 
################################################
#### Jingying Wang <astro.jywang@gmail.com> ####
################################################
