import numpy as np
import scipy.optimize as opt
import numpy.ma as ma
from . import timestream as kl
from .utils import legendre_timeseries, log_normal

Tcmb = 2.725
# for MeerKAT
d_freq = 208984.375
dump_period = 1.999154243


params = {
    'nd_ratio':         None,
    'ratio':            None,
    'Tptr':             None,
    'eta_p':            None,
    'Tnd':              None,
    'Tnd_ref':          None,
    'Tnd_std':          None,
    'Tel':              None,
    'Tgal':             None,
    'gain_coeffs':      None,
    'sm_coeffs':        None,
    'nd_0':             None,
    'nd_1a':            None,
    'nd_1b':            None
}


def gain_poly(t, gt):
    return legendre_timeseries(t, gt)


def func_sm(t, p_sm):
    return legendre_timeseries(t, p_sm)


# For track
def calc_total_model(timestamps, params):
    """
    Calculate total sky model as a function of time.

    Parameters
    ----------
    params : dict
        Dictionary of model parameters.

    Returns
    -------
    model : array_like
        Total model of the data as a function of time.
    """
    p = params
    gain_t = gain_poly(timestamps, p['gain_coeffs'])
    sm = func_sm(timestamps, p['sm_coeffs'])
    return gain_t * (sm + p['Tptr']*p['eta_p'] + p['Tel'] + p['Tgal']
                     + np.ones(len(timestamps))*p['Tcmb']
                     + vis_ndstamp_model(timestamps, p))


def calc_logprob(timestamps, vis, ch, params):
    """

    Parameters
    ----------

    Returns
    -------

    """
    p = params

    # Calculate total model and std. dev. estimate
    total_model = calc_total_model(timestamps, params)
    calc_error = total_model / gain_poly(timestamps, p['gain_coeffs']) \
        / np.sqrt(d_freq * dump_period)

    #
    result = ma.sum(log_normal(vis[:, ch], total_model, calc_error)) \
        + log_normal(p['Tnd'], p['Tnd_ref'], 0.1*p['Tnd_std']) \
        + log_normal(p['eta_p'], 1.0, 0.02)
    return result


def func_obj0(p, *args):
    """

    Parameters
    ----------

    Returns
    -------

    """
    timestamps, vis, ch, nd_ratio, Tptr, Tnd_ref, Tnd_std, Tel, Tgal, nd_0, nd_1a, nd_1b = args
    Tnd = p[0]
    eta_p = p[1]
    func_sm_param = p[2]
    func_gt_param = p[3:-1]
    ratio = p[-1]
    return -calc_logprob(timestamps, vis, ch, XXXX)


####
# def solve_params2(timestamps, vis_part1, vis_part2, ch, params):
    x0 = [Tnd_ref] + [eta_p0] + list(func_sm_param_a0) + list(func_sm_param_b0) \
        + list(func_gt_param0) + [ratio0]
    return opt.fmin_powell(func_obj2, x0=x0,
                           args=(timestamps, vis_part1, vis_part2, ch, nd_ratio, Tptr, Tnd_ref, Tnd_std, Tel, Tgal, nd_0, nd_1a, nd_1b))
####


def solve_params0(timestamps, vis, ch, nd_ratio, ratio0, Tptr, eta_p0, Tnd_ref, Tnd_std, Tel, Tgal, func_gt_param0, func_sm_param0, nd_0, nd_1a, nd_1b):
    """

    Parameters
    ----------

    Returns
    -------

    """
    x0 = [Tnd_ref]+[eta_p0]+list(func_sm_param0)+list(func_gt_param0)+[ratio0]
    return opt.fmin_powell(func_obj0,
                           x0=x0,
                           args=(timestamps, vis, ch, nd_ratio, Tptr, Tnd_ref, Tnd_std, Tel, Tgal, nd_0, nd_1a, nd_1b))


# -------------------------------------------------------------------------------
# For track, union fitting
# -------------------------------------------------------------------------------

def func_obj2(p, *args):
    """

    Parameters
    ----------

    Returns
    -------

    """
    timestamps, vis_part1, vis_part2, ch, nd_ratio, Tptr, Tnd_ref, Tnd_std, Tel, Tgal, nd_0, nd_1a, nd_1b = args
    Tnd = p[0]
    eta_p = p[1]
    func_sm_param_a = p[2]
    func_sm_param_b = p[3]
    func_gt_param = p[4:-1]
    ratio = p[-1:]
    delta_p = calc_logprob(timestamps, vis_part1, ch, XXX) \
        + calc_logprob(timestamps, vis_part2, ch, XXX)
    return -delta_p


def solve_params2(timestamps, vis_part1, vis_part2, ch, params):
    """

    Parameters
    ----------

    Returns
    -------

    """
    x0 = [Tnd_ref] + [eta_p0] + list(func_sm_param_a0) + list(func_sm_param_b0) \
        + list(func_gt_param0) + [ratio0]
    return opt.fmin_powell(func_obj2, x0=x0,
                           args=(timestamps, vis_part1, vis_part2, ch, nd_ratio, Tptr, Tnd_ref, Tnd_std, Tel, Tgal, nd_0, nd_1a, nd_1b))


# -------------------------------------------------------------------------------
# For scan
# -------------------------------------------------------------------------------

def solve_params_sm(timestamps, vis, ch, params):
    """

    Parameters
    ----------

    Returns
    -------

    """
    # Initial parameter values
    x0 = [eta_p0]+list(func_sm_param0)+list(func_gt_param0)+[ratio0]

    #####
    def func_obj_sm():
        timestamps, vis, ch, nd_ratio, Tptr, Tnd, Tel, Tgal, nd_0, nd_1a, nd_1b = args
        eta_p = p[0]
        func_sm_param = p[1:-6]
        func_gt_param = p[-6:-1]
        ratio = p[-1]

        return -calc_logprob_sm(timestamps, vis, ch, nd_ratio, ratio, Tptr, eta_p, Tnd, Tel, Tgal, func_gt_param, func_sm_param,
                                nd_0, nd_1a, nd_1b)
    #######

    return opt.fmin_powell(func_obj_sm,
                           x0=x0,
                           args=(timestamps, vis, ch, nd_ratio, Tptr, Tnd, Tel, Tgal,  nd_0, nd_1a, nd_1b))


def func_obj_sm(p, *args):
    """

    Parameters
    ----------

    Returns
    -------

    """


def calc_logprob_sm(timestamps, vis, ch, params):
    """

    Parameters
    ----------

    Returns
    -------

    """
    # Calculate total model and error
    total_model = calc_total_model(timestamps, params)
    calc_error = total_model / gain_poly(timestamps, p['gain_coeffs']) \
        / np.sqrt(d_freq*dump_period)

    # result = ma.sum(-(vis[:,ch]-total_model)**2/(2*error**2)-np.log(2*np.pi*error**2)/2.0) #supposing Gaussian

    # no point source at the moment
    result = ma.sum(log_normal(vis[:, ch], total_model, calc_error)) \
        + log_normal(params['eta_p'], 1.0, 1e-30)
    return result


def vis_ndstamp_model(timestamps, params):
    """

    Parameters
    ----------

    Returns
    -------

    """
    p = params
    vis_ndstamp = np.ma.array(np.zeros_like(timestamps))
    vis_ndstamp[nd_0] = 0
    vis_ndstamp[nd_1a] = Tnd * ratio
    vis_ndstamp[nd_1b] = Tnd * (nd_ratio - ratio)  # 1.8s injection
    return vis_ndstamp


def cal_gain0(fname, data, ant, pol, flags, ch,
              dp_tt, dp_ss, ang_deg, T_ptr, vis_clean):
    """

    Parameters
    ----------

    Returns
    -------

    """
    if fname in ['1551055211', '1551037708']:
        dp_ca, dp_cb, dp_c0a, dp_c1a, dp_c2a, dp_c3a, dp_c4a, dp_c0b, dp_c1b, dp_c2b, dp_c3b, dp_c4b = kl.cal_dp_c(fname, data, ant, pol, flags, ch, dp_tt, dp_ss, ang_deg)
        
        # vis gap for calibrator
        a1 = vis_clean[dp_c0a, ch].min()-vis_clean[dp_ca, ch].min()
        
        # T model gap for calibrator
        b1 = T_ptr[dp_ca].max() - T_ptr[dp_ca].min()
        
        # vis gap for calibrator
        a2 = vis_clean[dp_c0b, ch].min() - vis_clean[dp_cb, ch].min()
        
        # T model gap for calibrator
        b2 = T_ptr[dp_cb].max() - T_ptr[dp_cb].min()
    else:
        dp_ca, dp_cb, dp_c0a, dp_c1a, dp_c0b, dp_c1b = kl.cal_dp_c(
            fname, data, ant, pol, flags, ch, dp_tt, dp_ss, ang_deg)
        
        # vis gap for calibrator
        a1 = vis_clean[dp_c1a, ch].min() - vis_clean[dp_c0a, ch].min()
        
        # T model gap for calibrator
        b1 = T_ptr[dp_c1a].min() - T_ptr[dp_c0a].min()
        
        # vis gap for calibrator
        a2 = vis_clean[dp_c1b, ch].min() - vis_clean[dp_c0b, ch].min()
        
        # T model gap for calibrator
        b2 = T_ptr[dp_c1b].min() - T_ptr[dp_c0b].min()
    
    ga0 = a1 / b1
    gb0 = a2 / b2
    
    return ga0, gb0


