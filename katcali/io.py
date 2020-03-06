import katdal
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np
import pickle
from . import models as km
from .dataset import DataSet
import warnings

# Load default dataset
#ds = DataSet("meerkat2019.json")


#def load_data(fname):
#    warnings.warn("Deprecated; use DataSet object instead.", DeprecationWarning)
#    metadata = ds.get_metadata(fname)
#    return metadata


#def check_ants(fname):
#    warnings.warn("Deprecated; use DataSet object instead.", DeprecationWarning)
#    bad_ants = ds.get_bad_ants(fname)
#    calsrc = get_calibrator(fname)
#    target = calsrc.calname
#    cc = calsrc.coords
#    flux_model = calsrc.flux_model

#    print('calibrator: {target}, ra,dec= {ra}, {dec}'.format(
#        target=target, ra=cc.ra, dec=cc.dec))
#    print('bad_ants: %s' % str(bad_ants))
#    return target, cc, bad_ants, flux_model


#def ant_list(data):
#    warnings.warn("Deprecated; use DataSet object instead.", DeprecationWarning)
#    ant_list = []
#    for corr_i in range(len(data.ants)):
#        # check auto-corr
#        assert(data.corr_products[corr_i][0] == data.corr_products[corr_i][1])
#        ant_list.append(data.corr_products[corr_i][0][0:4])
#    return np.array(ant_list)


def cal_corr_id(data, recv):
    for i in range(len(data.ants)*2):
        a = data.corr_products[i]
        if (a[0] == recv and a[1] == recv):
            corr_id = i
            break
    return corr_id


#def load_ndparam(fname, data):
#    warnings.warn("Deprecated; use DataSet object instead.", DeprecationWarning)
#    nd_set = float(fname)
#    nd_time = 1.799235
#    nd_cycle = 19.9915424299  # only for stable diode noise pattern
#    nd_ratio = 1.8/data.dump_period  # 1.8???test
#    return nd_set, nd_time, nd_cycle, nd_ratio


#def load_tf(data):
#    warnings.warn("Deprecated; use DataSet object instead.", DeprecationWarning)
#    freqs = data.freqs
#    timestamps = data.timestamps
#    return timestamps, freqs


def load_ang_deg(ra, dec, cc):
    c_obs = SkyCoord(ra*u.deg,  dec*u.deg, frame='icrs')
    ang_deg = cc.separation(c_obs)/u.deg
    return ang_deg
