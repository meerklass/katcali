from seek.mitigation import sum_threshold
import numpy as np
import numpy.ma as ma
from . import filter as kf
from .timestream import time_indices #diode as kd
# param set ref: https://seek.readthedocs.io/en/latest/_modules/seek/mitigation/sum_threshold.html


def moving_average(x, n=3, axis=0):
    """
    Calculate a moving average of a 1D or 2D array using a uniform (boxcar) 
    window in one direction.
    
    This function returns an array the same shape as the input array, but there 
    are edge effects as a result (the window only has partial support near the 
    edges of the array).
    
    Parameters
    ----------
    x : array_like
        Input array.
    
    n : int, optional
        Size of boxcar window. Default: 3.
    
    axis : int, optional
        Which axis of the input array to apply the moving average to. 
        Default: 0.
    
    Returns
    -------
    avg : array_like
        Array with moving average applied. Has the same shape as input array.
    """
    # Make sure array is at least 2D
    if len(x.shape) == 1:
        x = x[:,np.newaxis]
    
    # Do moving average over chosen axis
    avg = np.zeros(x.shape, dtype=x.dtype)
    if axis == 1:
        for i in range(x.shape[0]):
            avg[i,:] = np.convolve(x[i,:], np.ones(n), mode='same') / float(n)
    else:
        for i in range(x.shape[1]):
            avg[:,i] = np.convolve(x[:,i], np.ones(n), mode='same') / float(n)
    return avg


def seek_rfi_mask(data, thres0, diode=False, verbose=False):
    """
    Use the "Hide and Seek" code (`sum_threshold`) to generate a set of RFI 
    flags.
    
    Parameters
    ----------
    data : array_like
        Time-ordered data, as a numpy masked array.
    
    thres0 : float
        First threshold to apply (chi_1).
    
    diode : bool, optional
        Whether noise diode is switched on in these data. If True, runs Seek 
        with di_kwargs = (1, 1); if False, runs Seek with di_kwargs = (25, 30).
    
    verbose : bool, optional
        Whether to print summary stats. Default: False.
    
    Returns
    -------
    masked_data : array_like
        Input data array with RFI flags added to mask.
    """
    # Settings for sum_threshold.get_rfi_mask
    sm_kwargs = sum_threshold.get_sm_kwargs(80, 80, 40, 40)
    if diode:
        # Noise diode is switched on
        di_kwargs = sum_threshold.get_di_kwrags(1, 1) # no expand
    else:
        # Noise diode is switched off
        di_kwargs = sum_threshold.get_di_kwrags(25, 30)
    
    # Get RFI mask
    rfi_mask = sum_threshold.get_rfi_mask( tod=data.astype('float'),
                                           mask=data.mask.astype('bool'),
                                           chi_1=thres0,
                                           sm_kwargs=sm_kwargs,
                                           di_kwargs=di_kwargs, 
                                           plotting=False )
    
    # Create numpy masked array with rfi_mask
    masked_data = ma.masked_array(data, mask=rfi_mask)
    
    # Print diagnostics
    if verbose:
        std_nomask = np.std(masked_data.compressed())
        std_nomask_norm = np.std(masked_data) / np.mean(masked_data)
        print('sigma(non-masked) = %2.2f' % std_nomask)
        print('sigma(non-masked) / mean(non-masked) = %2.2f' % std_nomask_norm)
    
    return masked_data


#def vis_flag(vis, flags, nd_label0, dp_w, First_Thresholds, verbose=False):
def vis_flag(vis, flags, idxs_nd_track, idxs_nd_scan, idxs_waste, thresholds, 
             flag_thres_time=0.5, flag_thres_freq=0.4, verbose=False):
    """
    Find and flag RFI in visibility data, using the Hide and Seek algorithm. 
    The noise diode on and off parts of the TOD are treated separately.
    
    Parameters
    ----------
    vis : array_like
        Time-ordered visibility data, as a 2D array.
    
    flags : array_like
        Array of existing/low-level flags on data (2D).
    
    idxs_nd_track : tuple of array_like
        Tuple of arrays of (raw, unflagged) time indices for samples taken in 
        tracking mode, split by whether noise diode is firing, in the order:
        
            track_on, track_off, track_on_start, track_on_stop
        
        Use `noise_diode_indices()` to generate these arrays in the appropriate 
        format
    
    idxs_nd_scan : tuple of array_like
        Tuple of arrays of (raw, unflagged) time indices for samples taken in 
        scan mode, split by whether noise diode is firing, in the order:
        
            scan_on, scan_off, scan_on_start, scan_on_stop
        
        Use `noise_diode_indices()` to generate these arrays in the appropriate 
        format.
    
    idxs_waste : array_like
        Array of time indices for wasted samples (use `timestream.time_indices`
        with label='waste' to get this array in the appropriate format).
    
    thresholds : tuple of float
        Initial thresholds (chi_1) to apply to each part of the TOD. This must 
        be a tuple of 4 floats, in the following order:
            
            thres_scan_off, thres_scan_on, thres_track_off, thres_track_on
        
        where 'on' and 'off' correspond to when the noise diode is on and off.
    
    flag_thres_time : float, optional
        If a time sample reaches a flag fraction above this value, flag all 
        frequency channels at this time. Default: 0.5.
    
    flag_thres_freq : float, optional
        If a frequency channel reaches a flag fraction above this value, flag 
        the whole channel for all times. Default: 0.4.
    
    verbose : bool, optional
        Whether to print summary stats and progress info. Default: False.
    
    Returns
    -------
    vis_final : array_like
        Input time-ordered visibility array with RFI flags added (returned as a 
        numpy masked array).
    """
    # Check that inputs have been passed in the right format
    assert isinstance(idxs_nd_track, tuple) and len(idxs_nd_track) == 4, \
        "'idxs_track' must be a tuple of 4 numpy arrays in the order: " \
        "track_on, track_off, track_on_start, track_on_stop"
    assert isinstance(idxs_nd_scan, tuple) and len(idxs_nd_scan) == 4, \
        "'idxs_scan' must be a tuple of 4 numpy arrays in the order: " \
        "scan_on, scan_off, scan_on_start, scan_on_stop"
    assert isinstance(idxs_waste, np.ndarray), \
        "'idxs_waste' must be a numpy array"
    assert isinstance(thresholds, tuple) and len(thresholds) == 4, \
        "'thresholds' must be a tuple of 4 floats in the order: " \
        "thres_scan_off, thres_scan_on, thres_track_off, thres_track_on"
    
    # Unpack tuples
    track_on, track_off, track_on_start, track_on_stop = idxs_nd_track
    scan_on, scan_off, scan_on_start, scan_on_stop = idxs_nd_scan
    thres_scan_off, thres_scan_on, thres_track_off, thres_track_on = thresholds
    
    # Create new masked array with low-level flags applied and mask waste data
    vis_masked = np.ma.array(vis, mask=flags)
    vis_masked[idxs_waste, :].mask = True
    
    # New masked array to add RFI mask to
    vis_clean = np.ma.array(np.zeros_like(vis_masked), mask=True)
    
    # Apply RFI flags to scan part of TOD
    if verbose: print("RFI flagging: scan TOD")
    vis_clean[scan_off,:]      = seek_rfi_mask(vis_masked[scan_off,:], 
                                               thres_scan_off, diode=False,
                                               verbose=verbose)
    vis_clean[scan_on_start,:] = seek_rfi_mask(vis_masked[scan_on_start,:], 
                                               thres_scan_on, diode=True,
                                               verbose=verbose)
    vis_clean[scan_on_stop,:]  = seek_rfi_mask(vis_masked[scan_on_stop,:], 
                                               thres_scan_on, diode=True, 
                                               verbose=verbose)
    
    # Apply RFI flags to track part of TOD
    if verbose: print("RFI flagging: track TOD")
    vis_clean[track_off,:]      = seek_rfi_mask(vis_masked[track_off,:], 
                                                thres_track_off, diode=False,
                                                verbose=verbose)
    vis_clean[track_on_start,:] = seek_rfi_mask(vis_masked[track_on_start,:], 
                                                thres_track_on, diode=True,
                                                verbose=verbose)
    vis_clean[track_on_stop,:]  = seek_rfi_mask(vis_masked[track_on_stop,:], 
                                                thres_track_on, diode=True,
                                                verbose=verbose)
    # Delete vis_masked to free-up memory
    del vis_masked
    
    # If the neighbouring diode-off sample is masked, the diode-on sample 
    # should be masked too
    for ch in range(np.shape(vis)[1]):
        for i in scan_on_start:
            if vis_clean.mask[i-1, ch]:
                vis_clean.mask[i, ch] = True

        for i in scan_on_stop:
            if vis_clean.mask[i+1, ch]:
                vis_clean.mask[i, ch] = True

        for i in track_on_start:
            if vis_clean.mask[i-1, ch]:
                vis_clean.mask[i, ch] = True

        for i in track_on_stop:
            if vis_clean.mask[i+1, ch]:
                vis_clean.mask[i, ch] = True
    
    # Completely flag frequencies and times with poor flag fractions
    vis_final = vis_clean.copy()
    t_len, ch_len = vis_clean.shape
    for i in range(ch_len):
        frac = float(np.array(vis_clean.mask[:,i] == True).sum()) / float(t_len)
        if frac > flag_thres_time:
            vis_final.mask[:,i] = True # flag all times at this freq
    for i in range(t_len):
        frac = float(np.array(vis_clean.mask[i,:] == True).sum()) / float(ch_len)
        if frac > flag_thres_freq:
            vis_final.mask[i,:] = True # flag all freqs at this time
    
    return vis_final


def mnad_flagger(vis, flags, thres=3., window=3, flag_thres_time=0.5, 
                 flag_thres_freq=0.4):
    """
    Median-normalised absolute difference (MNAD) filter, for flagging RFI. 
    
    This filter looks for RFI structure only in the frequency direction. It 
    first calculates the absolute value of the gradient of the data in the 
    frequency direction, then normalises this quantity by its median and 
    applies a threshold.
    
    The absolute value of the gradient can be smoothed using a boxcar moving 
    average, which helps pick out broader RFI features and reduce sensitivity 
    to noise fluctuations.
    
    Parameters
    ----------
    vis : array_like
        Array of visibility data, of shape (Ntimes, Nfreqs).
    
    flags : array_like
        Array of low-level flags, with the same shape as `vis`.
    
    thres : float, optional
        Threshold for MNAS statistic, above which a sample will be flagged. 
        Default: 3.
    
    window : int, optional
        Width of moving average boxcar in the frequency direction. Default: 3.
    
    flag_thres_time : float, optional
        If a time sample reaches a flag fraction above this value, flag all 
        frequency channels at this time. Default: 0.5.
    
    flag_thres_freq : float, optional
        If a frequency channel reaches a flag fraction above this value, flag 
        the whole channel for all times. Default: 0.4.
    
    Returns
    -------
    vis_masked : array_like
        Array of visibility data as a numpy masked array, with flags applied in 
        the mask.
    """
    # Apply low-level flags
    vis_masked = np.ma.masked_where(flags, vis) 
    
    # Calculate absolute difference in freq. direction
    abs_diff = np.abs(np.gradient(vis_masked, axis=1)) # abs. diff. in freq.
    
    # Calculate the median of the abs. diff. over freq. channels at each time
    median_per_time = np.zeros((vis_masked.shape[0],), dtype=vis_masked.dtype)
    for i in range(vis_masked.shape[0]):
        median_per_time[i] = np.median(abs_diff[i,:], axis=0)
    median_per_time[np.where(median_per_time == 0.)] = np.nan
    
    # Apply moving average to absolute differences (in freq. direction)
    smooth_abs_diff = moving_average(abs_diff, window, axis=1)
    
    # Calculate median-normalised abs. diff. statistic and apply mask
    mnad = smooth_abs_diff / np.atleast_2d(median_per_time).T
    vis_masked = np.ma.masked_where(mnad >= thres, vis_masked)
    
    # Mask frequencies and times with high flag occupancy
    frac_masked_freqs = np.mean(vis_masked.mask, axis=0)
    frac_masked_times = np.mean(vis_masked.mask, axis=1)
    for i in np.where(frac_masked_freqs > flag_thres_freq)[0]:
        vis_masked.mask[:,i] = True
    for i in np.where(frac_masked_times > flag_thres_time)[0]:
        vis_masked.mask[i,:] = True
    
    return vis_masked


############## single channel mask #######################################
def vis_flag_1ch(vis_clean, nd_labels, ch):
    nd_s1a, nd_s1b, nd_s1, nd_s0, nd_t1a, nd_t1b, nd_t1, nd_t0 = nd_labels
    print('group shape (with flags):')
    print(len(nd_s1a), len(nd_s1b), len(nd_s1), len(nd_s0),
          len(nd_t1a), len(nd_t1b), len(nd_t1), len(nd_t0))
    nd_s0_clean = []
    for i in nd_s0:
        if vis_clean.mask[i, ch] == False and vis_clean[i, ch] > 0:
            nd_s0_clean.append(i)

    for iter in range(5):
        print('iter='+str(iter))

        if iter == 0:
            nd_s0_i = nd_s0_clean
            l1 = len(nd_s0_clean)
        if iter > 0:
            nd_s0_i = nd_s0_clean2
            l1 = len(nd_s0_clean2)

        # output is the data want to keep
        nd_s0_clean2 = kf.curve_filter(nd_s0_i, vis_clean[nd_s0_i, ch])
        l2 = len(nd_s0_clean2)  # end len
        # print(l1,l2)

        if l1 == l2:
            break
    nd_s0_clean = nd_s0_clean2

    # update vis_clean
    for i in nd_s0:
        if i not in nd_s0_clean:
            #print i
            vis_clean.mask[i, ch] = True

    for i in nd_s1a:
        if vis_clean.mask[i-1, ch] == True:
            vis_clean.mask[i, ch] = True

    for i in nd_s1b:
        if vis_clean.mask[i+1, ch] == True:
            vis_clean.mask[i, ch] = True
    
    return vis_clean, nd_s0_clean
    
