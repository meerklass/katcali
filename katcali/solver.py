import numpy as np
from numpy.polynomial.legendre import legval
import scipy.optimize as opt
import numpy.ma as ma
from . import label_dump as kl

Tcmb=2.725
#for MeerKAT
d_freq=208984.375 
dump_period=1.999154243


def func_gt(t, gt):
    t_start=t[0]
    t_end=t[-1]
    x = (t - t_start) / (t_end - t_start) * 2 - 1  # time to [-1,1]
    return legval(x, gt)
    
def func_sm(t, p_sm):
    t_start=t[0]
    t_end=t[-1]
    x = (t - t_start) / (t_end - t_start) * 2 - 1  # time to [-1,1]
    return legval(x, p_sm)
    
def log_normal(x, mu, sigma):
    return -((x - mu)**2 /(2* sigma ** 2)) - np.log(np.sqrt(2 * np.pi) * np.abs(sigma))


#######################for track############################################
def calc_total_model(timestamps, nd_ratio, ratio, Tptr, eta_p, Tnd, Tel, Tgal, func_gt_param, func_sm_param, nd_0, nd_1a, nd_1b):
    return func_gt(timestamps,func_gt_param)*(func_sm(timestamps, func_sm_param)
                                              +Tptr*eta_p
                                              +Tel 
                                              +Tgal +np.ones(len(timestamps))*Tcmb
                                              +vis_ndstamp_model(timestamps, nd_0, nd_1a, nd_1b, Tnd, nd_ratio, ratio))


def calc_logprob(timestamps, vis, ch, nd_ratio, ratio, Tptr, eta_p, Tnd, Tnd_ref, Tnd_std, Tel, Tgal, func_gt_param, func_sm_param, nd_0, nd_1a, nd_1b):
    total_model=calc_total_model(timestamps, nd_ratio, ratio, Tptr, eta_p, Tnd, Tel, Tgal, func_gt_param, func_sm_param, nd_0, nd_1a, nd_1b)
    calc_error=total_model/func_gt(timestamps,func_gt_param)/np.sqrt(d_freq*dump_period)
   
    result=ma.sum(log_normal(vis[:,ch],total_model,calc_error))+log_normal(Tnd, Tnd_ref, 0.1*Tnd_std)+log_normal(eta_p, 1.0, 0.02)
    return result


def func_obj0(p, *args):
    timestamps, vis, ch, nd_ratio, Tptr, Tnd_ref, Tnd_std, Tel, Tgal, nd_0, nd_1a, nd_1b=args
    Tnd=p[0]
    eta_p=p[1]
    func_sm_param=p[2]
    func_gt_param=p[3:-1]
    ratio=p[-1]      
    return -calc_logprob(timestamps, vis, ch, nd_ratio, ratio, Tptr, eta_p, Tnd, Tnd_ref, Tnd_std, Tel, Tgal, func_gt_param, func_sm_param, nd_0, nd_1a, nd_1b)


def solve_params0(timestamps, vis, ch, nd_ratio, ratio0, Tptr, eta_p0, Tnd_ref, Tnd_std, Tel, Tgal, func_gt_param0, func_sm_param0, nd_0, nd_1a, nd_1b):
    return opt.fmin_powell(func_obj0, 
                           x0=[Tnd_ref]+[eta_p0]+list(func_sm_param0)+list(func_gt_param0)+[ratio0],
                           args=(timestamps, vis, ch, nd_ratio, Tptr, Tnd_ref, Tnd_std, Tel, Tgal, nd_0, nd_1a, nd_1b))

####for track, union fitting########
def func_obj2(p, *args):
    timestamps, vis_part1, vis_part2, ch, nd_ratio, Tptr, Tnd_ref, Tnd_std, Tel, Tgal, nd_0, nd_1a, nd_1b=args
    Tnd=p[0]
    eta_p=p[1]
    func_sm_param_a=p[2]
    func_sm_param_b=p[3]
    func_gt_param=p[4:-1]
    ratio=p[-1:]      
    return -calc_logprob(timestamps, vis_part1, ch, nd_ratio, ratio, Tptr, eta_p, Tnd, Tnd_ref, Tnd_std, Tel, Tgal, func_gt_param, func_sm_param_a, nd_0, nd_1a, nd_1b)-calc_logprob(timestamps, vis_part2, ch, nd_ratio, ratio, Tptr, eta_p, Tnd,  Tnd_ref, Tnd_std, Tel, Tgal, func_gt_param, func_sm_param_b, nd_0, nd_1a, nd_1b)


def solve_params2(timestamps, vis_part1, vis_part2, ch, nd_ratio, ratio0, Tptr, eta_p0, Tnd_ref,Tnd_std, Tel, Tgal, func_gt_param0, func_sm_param_a0,  func_sm_param_b0, nd_0, nd_1a, nd_1b):
    return opt.fmin_powell(func_obj2, 
                           x0=[Tnd_ref]+[eta_p0]+list(func_sm_param_a0)+list(func_sm_param_b0)+list(func_gt_param0)+[ratio0],
                           args=(timestamps, vis_part1, vis_part2, ch, nd_ratio, Tptr, Tnd_ref, Tnd_std, Tel, Tgal, nd_0, nd_1a, nd_1b))



#################for scan#############################
def solve_params_sm(timestamps, vis, ch, nd_ratio, ratio0, Tptr, eta_p0, Tnd, Tel, Tgal, func_gt_param0, func_sm_param0,
                      nd_0, nd_1a, nd_1b):
    return opt.fmin_powell(func_obj_sm, 
                           x0=[eta_p0]+list(func_sm_param0)+list(func_gt_param0)+[ratio0], 
                           args=(timestamps, vis, ch, nd_ratio, Tptr, Tnd, Tel, Tgal,  nd_0, nd_1a, nd_1b))

def func_obj_sm(p, *args):
    timestamps, vis, ch, nd_ratio, Tptr, Tnd, Tel, Tgal, nd_0, nd_1a, nd_1b=args
    eta_p=p[0]
    func_sm_param=p[1:-6]
    func_gt_param=p[-6:-1]
    ratio=p[-1]
    
    return -calc_logprob_sm(timestamps, vis, ch, nd_ratio, ratio, Tptr, eta_p, Tnd, Tel, Tgal, func_gt_param, func_sm_param,
                              nd_0, nd_1a, nd_1b)


def calc_logprob_sm(timestamps, vis, ch, nd_ratio, ratio, Tptr, eta_p, Tnd, Tel,Tgal, func_gt_param, func_sm_param, nd_0, nd_1a, nd_1b):
    total_model=calc_total_model_sm(timestamps, nd_ratio, ratio, Tptr, eta_p, Tnd, Tel, Tgal, func_gt_param, func_sm_param, nd_0, nd_1a, nd_1b)
    calc_error=total_model/func_gt(timestamps,func_gt_param)/np.sqrt(d_freq*dump_period)
    #result = ma.sum(-(vis[:,ch]-total_model)**2/(2*error**2)-np.log(2*np.pi*error**2)/2.0) #supposing Gaussian
    result=ma.sum(log_normal(vis[:,ch],total_model,calc_error))+log_normal(eta_p, 1.0, 1e-30) #no point source at the moment
    return result

def calc_total_model_sm(timestamps, nd_ratio, ratio, Tptr, eta_p, Tnd, Tel, Tgal, func_gt_param, func_sm_param, nd_0, nd_1a, nd_1b):
    
    return func_gt(timestamps,func_gt_param)*(func_sm(timestamps, func_sm_param)+Tel
                                              +Tgal+np.ones(len(timestamps))*Tcmb ###Galactic and CMB added!
                                              +Tptr*eta_p
                                              +vis_ndstamp_model(timestamps, nd_0, nd_1a, nd_1b, Tnd, nd_ratio, ratio))




def vis_ndstamp_model(timestamps,nd_0,nd_1a,nd_1b, Tnd, nd_ratio, ratio):
    vis_ndstamp=np.ma.array(np.zeros_like(timestamps))
    vis_ndstamp[nd_0]=0
    vis_ndstamp[nd_1a]=Tnd*ratio
    vis_ndstamp[nd_1b]=Tnd*(nd_ratio-ratio) #1.8s injection
    return vis_ndstamp

def cal_gain0(fname,data,ant,pol,flags,ch,dp_tt,dp_ss,ang_deg,T_ptr,vis_clean):
    if fname in ['1551055211','1551037708']:
        dp_ca,dp_cb,dp_c0a, dp_c1a,dp_c2a,dp_c3a,dp_c4a,dp_c0b,dp_c1b,dp_c2b,dp_c3b,dp_c4b=kl.cal_dp_c(fname,data,ant,pol,flags,ch,dp_tt,dp_ss,ang_deg)
        a1=vis_clean[dp_c0a,ch].min()-vis_clean[dp_ca,ch].min() #vis gap for calibrator
        b1=T_ptr[dp_ca].max()-T_ptr[dp_ca].min() # T model gap for calibrator
        a2=vis_clean[dp_c0b,ch].min()-vis_clean[dp_cb,ch].min() #vis gap for calibrator
        b2=T_ptr[dp_cb].max()-T_ptr[dp_cb].min() # T model gap for calibrator
    else:
        dp_ca,dp_cb,dp_c0a, dp_c1a,dp_c0b,dp_c1b=kl.cal_dp_c(fname,data,ant,pol,flags,ch,dp_tt,dp_ss,ang_deg)
        a1=vis_clean[dp_c1a,ch].min()-vis_clean[dp_c0a,ch].min() #vis gap for calibrator
        b1=T_ptr[dp_c1a].min()-T_ptr[dp_c0a].min() # T model gap for calibrator
        a2=vis_clean[dp_c1b,ch].min()-vis_clean[dp_c0b,ch].min() #vis gap for calibrator
        b2=T_ptr[dp_c1b].min()-T_ptr[dp_c0b].min() # T model gap for calibrator
    ga0=a1/b1
    gb0=a2/b2
    return ga0,gb0




'''
#######only for test
########################for scan, break bkg for each line##############

def solve_params_smbr(timestamps, vis, ch, ratio0, Tptr, Tnd, func_gt_param0, func_sm_param0,
                      nd_0, nd_1a, nd_1b, nt_az_edge, error):
    return opt.fmin_powell(func_obj_smbr, 
                           x0=list(func_sm_param0)+list(func_gt_param0)+[ratio0], 
                           args=(timestamps, vis, ch, Tptr, Tnd, nd_0, nd_1a, nd_1b, nt_az_edge, error))

def func_obj_smbr(p, *args):
    timestamps, vis, ch, Tptr, Tnd, nd_0, nd_1a, nd_1b, nt_az_edge, error=args
    func_sm_param=p[0:-5]
    func_gt_param=p[-5:-1]
    ratio=p[-1]
    
    return -calc_logprob_smbr(timestamps, vis, ch, ratio, Tptr, Tnd, func_gt_param, func_sm_param,
                              nd_0, nd_1a, nd_1b, nt_az_edge, error)



def calc_logprob_smbr(timestamps, vis, ch, ratio, Tptr, Tnd, func_gt_param, func_sm_param, nd_0, nd_1a, nd_1b, nt_az_edge, error):
    total_model=calc_total_model_smbr(timestamps, ratio, Tptr, Tnd, func_gt_param, func_sm_param, nd_0, nd_1a, nd_1b,nt_az_edge)
    result = ma.sum(-(vis[:,ch]-total_model)**2/(2*error**2)-np.log(2*np.pi*error**2)/2.0) 
    return result

def calc_total_model_smbr(timestamps, ratio, Tptr, Tnd, func_gt_param, func_sm_param, nd_0, nd_1a, nd_1b,nt_az_edge):
    
    return func_gt(timestamps,func_gt_param)*(func_sm_break_long(timestamps, func_sm_param, nt_az_edge)
                                              +gal[:,ch]+np.ones(len(timestamps))*Tcmb. ###Galactic and CMB added!
                                              +Tptr
                                              +vis_ndstamp_model(timestamps, nd_0, nd_1a, nd_1b, Tnd, ratio))

def func_sm_break_long(timestamps, func_sm_param, nt_az_edge):
    sm_scan=func_sm_break(timestamps, func_sm_param, nt_az_edge)
    return np.concatenate((np.ones(nt_az_edge[0])*sm_scan[0],
                                                sm_scan,
                                                np.ones(len(timestamps)-nt_az_edge[-1])*sm_scan[-1])) 

def func_sm_break(ts, p_sm, nt_az_edge): #ts is full list of timestamps
    from numpy.polynomial.legendre import legval
    sm_part=[]
    
    n=int(len(p_sm)/(len(nt_az_edge)-1))
    for i in range(len(nt_az_edge)-1):
        t_start=ts[nt_az_edge[i]]
        t_end=ts[nt_az_edge[i+1]] #not included
        x = (ts[nt_az_edge[i]:nt_az_edge[i+1]] - t_start) / (t_end - t_start) * 2 - 1  # time to [-1,1]
        sm_part+=list(legval(x, p_sm[n*i:n*(i+1)]))
    
    return np.array(sm_part)
'''