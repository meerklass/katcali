import numpy as np
from scipy import stats
from . import filter as kfilter
from .utils import valid_filename


def select_track(data, ant, pol):
    """

    Parameters
    ----------
    data : metadata object
        x

    ant : str
        Antenna, e.g. 'm006'.

    pol : str
        Polarisation, e.g. 'h'.

    Returns
    -------
    idxs_track : ?
        ?
    """
    data.select(ants=ant, pol=pol, scans='track')
    scans_t = []
    scans_tw = []
    for s in data.scans():
        if data.shape[0] > 50: # only need 2 mins
            scans_t.append(data.scan_indices[0])
        else:
            scans_tw.append(data.scan_indices[0])
    data.select(ants=ant, pol=pol, scans=scans_t)
    dp_t = data.dumps
    data.select()  # recover after select!!!
    data.select(ants=ant, pol=pol)
    return dp_t, scans_tw


def select_scan(data, ant, pol):
    """

    Parameters
    ----------
    data : metadata object
        x

    ant : str
        Antenna, e.g. 'm006'.

    pol : str
        Polarisation, e.g. 'h'.

    Returns
    -------
    idxs_waste : ?
        ?
    """
    data.select(ants=ant, pol=pol, scans='scan')
    scans_s = []
    scans_sw = []
    for s in data.scans():
        if data.shape[0] > 50:
            scans_s.append(data.scan_indices[0])
        else:
            scans_sw.append(data.scan_indices[0])
    data.select(ants=ant, pol=pol, scans=scans_s)
    dp_s = data.dumps
    # dp_sb=data.dumps[0]
    # dp_se=data.dumps[-1]
    data.select()  # recover after select!!!
    data.select(ants=ant, pol=pol)
    return dp_s, scans_sw


def select_waste(data, ant, pol):
    """

    Parameters
    ----------
    data : metadata object
        x

    ant : str
        Antenna, e.g. 'm006'.

    pol : str
        Polarisation, e.g. 'h'.

    Returns
    -------
    idxs_waste : ?
        ?
    """
    data.select(ants=ant, pol=pol, scans=('slew', 'stop'))
    dp_w1 = data.dumps

    dp_t, scans_tw = select_track(data, ant, pol)
    dp_s, scans_sw = select_scan(data, ant, pol)

    data.select(ants=ant, pol=pol, scans=scans_tw+scans_sw)
    dp_w2 = data.dumps
    dp_w = list(dp_w1)+list(dp_w2)
    dp_w.sort()
    dp_w = np.array(dp_w)
    data.select()  # recover after select!!!
    data.select(ants=ant, pol=pol)
    return dp_w


def time_indices(meta, ant, pol, label, flags=None):
    """
    Return an array of indices of time samples with a given label, e.g. 
    only time samples that were taken while the telescope was scanning.
    
    Parameters
    ----------
    meta : obj
        Metadata object returned by DataSet.get_metadata().
    
    ant : str
        Name of antenna, e.g. 'm006'.
    
    pol : str
        Name of polarisation, e.g. 'h' or 'v'.
    
    label : str
        Which label to find indices for. Options are:
         - track:       Telescope is tracking a source (flags applied)
         - scan:        Telescope is scanning (flags applied)
         - flagged:     Time sample was flagged
         - waste:       Wasted due to poor conditions
         - track-raw:   Telescope is tracking a source (flags ignored)
         - scan-raw:    Telescope is scanning (flags ignored)
    
    flags : array_like, optional
        2D array (waterfall) of flags. Only needed if label is not a 'raw' mode.
    
    Returns
    -------
    indices : array_like
        Indices of time samples that match `label` in the time array.
    """
    # Check that label is valid
    valid_labels = ['track', 'scan', 'flagged', 'waste', 'track-raw', 
                    'scan-raw']
    if label not in valid_labels:
        raise ValueError("'%s' is not a valid label type." % label)
    
    # Load flags if needed
    if 'raw' not in label and 'waste' not in label:
        assert flags is not None, "'flags' must be specified."
    
    # Return time indices with a given label
    if 'track' in label:
        idxs_track, scans_tw = select_track(meta, ant, pol)
        
        if 'raw' in label:
            # Leave flagged samples in data
            return idxs_track
        else:
            # Remove flagged samples
            f = flags[idxs_track]
            return idxs_track[~f] # negation of flags
            
    elif 'scan' in label:
        idxs_scan, scans_sw = select_scan(meta, ant, pol)
        
        if 'raw' in label:
            # Leave flagged samples in data
            return idxs_scan
        else:
            # Remove flagged samples
            f = flags[idxs_scan]
            return idxs_scan[~f] # negation of flags
            
    elif label == 'waste':
        idxs_waste = select_waste(meta, ant, pol)
        return idxs_waste
        
    elif label == 'flagged':
        return np.where(flags)[0]
        
    else:
        return # should never get here


def noise_diode_indices(times, offset, period, intersect=None):
    """
    Return the indices of time samples when the noise diode is ON and OFF.
    
    Parameters
    ----------
    times : array_like
        Array of time samples. Assumed to be regularly spaced.
    
    offset : int
        Index of first noise diode fire in the time series.
    
    period : int
        Number of time samples in between noise diode fires.
    
    intersect : array_like, optional
        If specified, return the intersections of this array with the output 
        arrays. Should be an integer array of time indices. Default: None.
    
    Returns
    -------
    idxs_on : array_like
        Time indices when noise diode is on.
    
    idxs_off : array_like
        Time indices when noise diode is off.
    
    idxs_on_start : array_like
        Time indices when noise diode first switches on. (Noise 
        diode fires span two time samples; this is the first in 
        each firing.)
    
    idxs_on_stop : array_like
        Time indices when noise diode stops being on. (Noise 
        diode fires span two time samples; this is the last in 
        each firing.)
    """
    assert len(times.shape) == 1, "'times' must be a 1D array"
    idxs = np.arange(times.size, dtype=np.integer)
    
    # Time indices when noise diode starts being on and ends being on
    # (noise diode fires span two time samples; the first is 'start' 
    # and second is 'stop')
    idxs_on_start = idxs[offset::period]
    idxs_on_stop = idxs[offset+1::period]
    
    # Indices of all time samples where the noise diode is on (union of 
    # on_start and on_end); interleave so that ordering should be correct
    idxs_on = np.zeros(idxs_on_start.size + idxs_on_stop.size, dtype=np.integer)
    idxs_on[::2] = idxs_on_start
    idxs_on[1::2] = idxs_on_stop
    
    # Time indices where noise diode is off
    idxs_off = idxs[~np.isin(idxs, idxs_on)]
    
    # Perform intersections
    if intersect is not None:
        idxs_on = np.intersect1d(idxs_on, intersect)
        idxs_off = np.intersect1d(idxs_off, intersect)
        idxs_on_start = np.intersect1d(idxs_on_start, intersect)
        idxs_on_stop = np.intersect1d(idxs_on_stop, intersect)
    
    return idxs_on, idxs_off, idxs_on_start, idxs_on_stop


def calibrator_time_indices(meta, ant, pol, idxs_track, idxs_scan, flags, ang_deg, 
                            targets=[0,1], sigma=10., niter=3):
    """
    Return the indices of time samples at various pointings on and around the 
    calibration source (pointing strategy varies for different observations).
    
    Parameters
    ----------
    meta : obj
        katdal metadata object.
    
    ant : str
        Antenna to use.
        
    pol : str
        Polarisation to use.
    
    idxs_track, idxs_scan : array_like
        Time indices for samples taken during tracking and scanning modes 
        respectively.
        
        This function finds the calibration source time indices from the 
        tracking time indices. The scanning indices are only needed to split 
        the cal. source time indices into pre-scanning and post-scanning groups.
    
    flags : array_like
        1D boolean array of flags for time samples. Should be for a single 
        frequency channel.
    
    ang_deg : float
        ???.
    
    targets : list
        List of target field IDs to cycle through. Not every file has the same 
        calibration strategy; all files have targets 0 (on-source) and 1 (off-
        source, in outskirts), while some have several other off-source 
        pointings.
    
    sigma : float, optional
        Default: 10.
    
    niter : int, optional
        Default: 3.
    
    Returns
    -------
    indices : dict
        Dictionary of time indices corresponding to calibration source 
        pointings.
        
        The keys of the dict correspond to the numerical ID of each calibration 
        pointing, e.g. 0.
        
        The items of the dict are tuples containing a pair of arrays, 
        (idxs_prescan, idxs_postscan), for calibration observations taken 
        before and after the scanning phase.
    """
    # Get start/stop time indices of scan
    idx_scan_start = idxs_scan[0]
    idx_scan_stop = idxs_scan[-1]
    
    # Loop over targets and get corresponding time indices
    calibrator_time_idxs = {}
    for target in targets:
        
        # Select target from tracking-mode observations
        meta.select()
        meta.select(ants=ant, pol=pol, scans='track', targets=target)
        
        # Select time samples that were taken in track mode and aren't flagged
        flag_segment = flags[meta.dumps] # get flags for this timestream segment
        idxs_cal = np.intersect1d(meta.dumps[~flag_segment], idxs_track)
        
        idxs_cal = kfilter.deg_filter(idxs_cal, ang_deg=ang_deg, sigma=sigma, 
                                      niter=niter)
        idxs_prescan = idxs_cal[idxs_cal < idx_scan_start] # before scan starts
        idxs_postscan = idxs_cal[idxs_cal > idx_scan_stop] # after scan stops
        
        # Reset select (FIXME: assumed to not be a time-consuming operation)
        meta.select()
        
        # Add to dict
        idx_dict[target] = (idxs_prescan, idxs_postscan)
        
    #if fname in ['1551055211', '1551037708'], do 2,3,4
    return idx_dict
    

