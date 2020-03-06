import numpy as np
import json
import pickle
import katdal
from .models import cal_sources
from .utils import valid_filename, file_size
from . import timestream


class DataSet(object):
    
    def __init__(self, dbname):
        """
        Contains a list of datafiles to use, along with their metadata.
        
        Each file 
        
        Parameters
        ----------
        dbname : str
            File containing list of datasets.
        """
        self.dbname = dbname
        self.files = {} # info about file
        self._metadata = {}
        try:
            self.load()
        except:
            raise
    
    
    def load(self):
        """
        Load info about dataset from file.
        """
        with open(self.dbname, 'r') as f:
            self.files = json.load(f)
    
    
    def save(self, dbname=None):
        """
        Save info about dataset to a file.
        
        Parameters
        ----------
        dbname : str, optional
            If specified, save dataset info to this filename. 
            Default: None (save to original file).
        """
        if dbname is None: dbname = self.dbname
        with open(dbname, 'w') as outfile:
            json.dump(self.files, outfile, indent=2)
    
    
    def load_metadata(self, filename, verbose=False):
        """
        Load metadata into memory for a given filename.
        
        Parameters
        ----------
        filename : str
            Name of metadata file to load (in JSON format).
        
        verbose : bool, optional
            Whether to print simple progress messages during load. 
            Default: True.
        """
        # Validate filename
        filename = valid_filename(filename)
        
        # Check if data are already loaded for this file
        if filename in self._metadata.keys():
            if verbose:
                print("Metadata already loaded for '%s'." % filename)
            return
        
        # Get info about this file
        info = self.files[filename]
        
        # Construct metadata filename
        tmpl = "{proj_root}/{data_root}/{fname}/{fname}/{fname}_sdp_l0.full.rdb"
        fname = tmpl.format(proj_root=info['project_root'], 
                            data_root=info['data_root'],
                            fname=filename)
        
        # Check file size
        fsize, units = file_size(fname)
        if verbose: print("Loading metadata file (%3.1f %s)" % (fsize, units))
        
        # Load the metadata
        meta = katdal.open(fname)
        self._metadata[filename] = meta
    
    
    def load_data(self, filename, ant, pol, verbose=True):
        return self.get_data(filename, ant, pol, verbose=verbose)
    
    
    def get_data(self, filename, ant, pol, verbose=True):
        """
        Get visibility data from a given file for a given receiver.
        
        Parameters
        ----------
        filename : str
            Name of file to load data for.
        
        ant : str
            Name of antenna to load data for, e.g. 'm006'.
        
        pol : str
            Name of polarisation to load data for, e.g. 'h' or 'v'.
        
        verbose : bool, optional
            Print simple status messages while loading data. Default: True.
        
        Returns
        -------
        vis : ndarray
            Visibility data.
        
        flags : ndarray
            Flags for visibility data.
        """
        # Validate filename
        filename = valid_filename(filename)
        
        # Get info about this file
        info = self.files[filename]
        
        # Construct receiver-polarisation string
        recvp = "%s%s" % (ant, pol)
        
        # Construct filename
        tmpl = "{proj_root}/raw_vis/{data_root}/{fname}/{fname}_{recvp}_vis_data"
        fname = tmpl.format(proj_root=info['project_root'], 
                            data_root=info['data_root'],
                            fname=filename,
                            recvp=recvp)
        
        # Check file size
        fsize, units = file_size(fname)
        print("Data file: %3.1f %s" % (fsize, units))
        
        # Load data using pickle
        with open(fname, 'rb') as f:
            # FIXME: Encoding needed for Py2 pickles
            data = pickle.load(f, encoding='latin1')
            
            # Test that correct receiver has been found
            recv_pair = data['recv_pair']
            rp0 = recv_pair[0]
            if isinstance(rp0, bytes):
                # Handle bytes vs string issues between Py2 and Py3
                rp0 = rp0.decode()
            assert(rp0 == recvp), \
                "Unexpected receiver/pol. '%s' found in file '%s'" \
                % (rp0, fname)
            
            # Get visibility and flag data
            vis = data['vis']
            flags = data['flags']
            
        return vis, flags
    
    
    def get_metadata(self, filename):
        """
        Return katdal metadata object for a given file.
        
        Parameters
        ----------
        filename : str
            Name of file.
        
        Returns
        -------
        metadata : katdal object
            katdal metadata object.
        """
        # Validate filename
        filename = valid_filename(filename)
        
        # Make sure metadata is loaded
        self.load_metadata(filename)
        return self._metadata[filename]
    
    
    def get_ants(self, filename):
        """
        Return a list of available antennas (which have data) for a given file.
        
        N.B. This method will load the metadata for a file if it is not already 
        in memory. This can take a few seconds.
        
        Parameters
        ----------
        filename : str
            Name of file.
        
        Returns
        -------
        ants : list of str
            List of antennas.
        """
        # Validate filename
        filename = valid_filename(filename)
        
        # Make sure metadata is loaded
        meta = self.get_metadata(filename)
        
        # Loop over available antennas
        ant_list = []
        for corr_i in range(len(meta.ants)):
            prods = meta.corr_products[corr_i]
            assert prods[0] == prods[1] # check that this is auto-corr only
            ant_list.append(prods[0][0:4])
        return ant_list
    
    
    def get_bad_ants(self, filename):
        """
        Return a list of bad antennas for a given file.
        
        Parameters
        ----------
        filename : str
            Name of file.
        
        Returns
        -------
        bad_ants : list of str
            List of bad antennas.
        """
        # Validate filename
        filename = valid_filename(filename)
        
        # Return list of bad antennas
        return self.files[filename]['bad_ants']
    
    
    def get_calibrator(self, filename):
        """
        Return the source model CalSource object corresponding to the 
        calibration source that was used for a given file.
        
        Parameters
        ----------
        filename : str
            Name of file to find calibrator source model for.
        
        Returns
        -------
        cal : CalSource
            Returns instance of CalSource class corresponding to the primary 
            calibration source for this file.
        """
        # Validate filename
        filename = valid_filename(filename)
        
        # Get calibrator for this field
        target = self.files[filename]['target']
        if target in cal_sources.keys():
            return cal_sources[target]
        else:
            raise KeyError("File '%s' has calibrator source '%s', but a model "
                           "for this source was not found." % (filename, target))
    
    
    def get_time_indices(self, filename, ant, pol, label, ch=None):
        """
        Return an array of indices of time samples with a given label, e.g. 
        only time samples that were taken while the telescope was scanning.
        
        Parameters
        ----------
        filename : str
            Name of file.
        
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
        
        ch : int, optional
            Which channel to retrieve flags for, if flagging is to be applied. 
            Default: None.
        
        Returns
        -------
        indices : array_like
            Indices of time samples that match `label` in the time array.
        """
        # Validate filename
        filename = valid_filename(filename)
        
        # Make sure metadata is loaded
        meta = self.get_metadata(filename)
        
        # Load flags if needed
        if 'raw' not in label:
            assert ch is not None, "Channel 'ch' must be specified."
            _vis, flags = self.get_data(filename=filename, ant=ant, pol=pol)
            flags = flags[:,ch]
        
        # Get time indices for this label
        return timestream.time_indices(meta=meta, ant=ant, pol=pol, 
                                       label=label, flags=flags)
    
    
    def get_calibrator_indices(self, filename, ant, pol, ch, ang_deg, 
                               targets=[0,1], sigma=10., niter=3):
        """
        Return an array of indices of time samples during observations of the 
        calibration source.
        
        Parameters
        ----------
        filename : str
            Name of file.
        
        ant : str
            Name of antenna, e.g. 'm006'.
        
        pol : str
            Name of polarisation, e.g. 'h' or 'v'.
        
        ch : int
            Which channel to retrieve flags for.
        
        ang_deg : float
            ???.
        
        targets : list
            List of target field IDs to cycle through. Not every file has the 
            same calibration strategy; all files have targets 0 (on-source) 
            and 1 (off-source, in outskirts), while some have several other 
            off-source pointings.
        
        sigma : float, optional
            Default: 10.
        
        niter : int, optional
            Default: 3.
        
        Returns
        -------
        indices : dict
            Dictionary of time indices corresponding to calibration source 
            pointings.
            
            The keys of the dict correspond to the numerical ID of each 
            calibration pointing, e.g. 0.
            
            The items of the dict are tuples containing a pair of arrays, 
            (idxs_prescan, idxs_postscan), for calibration observations taken 
            before and after the scanning phase.
        """
        # Validate filename
        filename = valid_filename(filename)
        
        # Make sure metadata is loaded
        meta = self.get_metadata(filename)
        
        # Load flags
        _vis, flags = self.get_data(filename=filename, ant=ant, pol=pol)
        flags = flags[:,ch]
        
        # Get time indices for tracking and scanning modes
        idxs_track = timestream.time_indices(meta=meta, ant=ant, pol=pol, 
                                             label='track-raw', flags=flags)
        idxs_scan = timestream.time_indices(meta=meta, ant=ant, pol=pol, 
                                            label='scan-raw', flags=flags)
        
        # Get time indices and return as a dict
        return timestream.calibrator_time_indices(meta, idxs_track, idxs_scan, 
                                                  flags=flags, ang_deg=ang_deg, 
                                                  targets=targets, sigma=sigma, 
                                                  niter=niter)    
    
    
    def get_noise_diode_indices(self, filename, intersect):
        """
        Return an array of indices of time samples with a given label, e.g. 
        only time samples that were taken while the telescope was scanning.
        
        Parameters
        ----------
        filename : str
            Name of file.
            
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
        times = self.get_times(filename)
        
        # Get noise diode settings
        offset = self.get_diode_info(filename, field='offset')
        period = self.get_diode_info(filename, field='jump_limit')
        
        return timestream.noise_diode_indices(times=times, offset=offset, 
                                              period=period, 
                                              intersect=intersect)
    
    
    def get_coords(self, filename):
        """
        Return arrays of RA, Dec, azimuth, and elevation for each time sample.
        
        N.B. This method will load the metadata for a file if it is not already 
        in memory. This can take a few seconds.
        
        Parameters
        ----------
        filename : str
            Name of file.
        
        Returns
        -------
        ra, dec, az, el : array_like
            Arrays of RA, Dec, Az, and El values for each time sample.
        """
        # Validate filename
        filename = valid_filename(filename)
        
        # Make sure metadata is loaded
        meta = self.get_metadata(filename)
        
        # Get arrays
        ra = meta.ra[:,0]
        dec = meta.dec[:,0]    
        az = meta.az[:,0]
        el = meta.el[:,0]
        
        return ra, dec, az, el
    
    
    def get_freqs(self, filename):
        """
        Return an array of frequency values for a given file.
        
        N.B. This method will load the metadata for a file if it is not already 
        in memory. This can take a few seconds.
        
        Parameters
        ----------
        filename : str
            Name of file to return frequency array for.
        
        Returns
        -------
        freqs : array_like
            Array of frequencies.
        """
        # Validate filename
        filename = valid_filename(filename)
        
        # Make sure metadata is loaded
        meta = self.get_metadata(filename)
        return meta.freqs
    
    
    def get_times(self, filename):
        """
        Return an array of timestamps for a given file.
        
        N.B. This method will load the metadata for a file if it is not already 
        in memory. This can take a few seconds.
        
        Parameters
        ----------
        filename : str
            Name of file to return timestamp array for.
        
        Returns
        -------
        times : array_like
            Array of timestamps.
        """
        # Validate filename
        filename = valid_filename(filename)
        
        # Make sure metadata is loaded
        meta = self.get_metadata(filename)
        return meta.timestamps
    
    
    def get_diode_info(self, filename, field=None):
        """
        Get information about the noise diode settings for a particular file.
        
        Several special fields are available that do not appear in the returned 
        `diode_info` dict, but which can be obtained by specifying the 
        following field names: nd_set, nd_time, nd_cycle, nd_ratio.
        
        Parameters
        ----------
        filename : str
            Name of file to find noise diode model for.
        
        field : str, optional
            If specified, return the data for a given setting only. Otherwise, 
            a dictionary of all settings for the file will be returned 
            (currently excluding special fields; see above).
            Default: None (return full dict).
        
        Returns
        -------
        diode_info : dict or dict item
            Dict of all noise diode settings for this file (if 'field' is not 
            specified). Otherwise, returns the dict item corresponding to 
            'field'.
        """
        # Validate filename
        filename = valid_filename(filename)
        
        # Special fields
        if field == 'nd_set':
            return float(filename)
        elif field == 'nd_time':
            return 1.799235
        elif field == 'nd_cycle':
            return 19.9915424299 # only for stable diode noise pattern
        elif field == 'nd_ratio':
            meta = self.get_metadata(filename)
            return 1.8 / meta.dump_period
        else:
            pass
        
        # Return whole dict if field is None
        if field is None:
            return self.files[filename]['diode']
        else:
            return self.files[filename]['diode'][field]
        
    
    def add_file(self, filename, target=None, root=None, diode=None, 
                 bad_ants=[]):
        """
        Add a file to the list of files within this dataset.
        
        Parameters
        ----------
        target : str
            Name of the calibration source target that was used for this file.
        
        root : str
            Root directory for the data and metadata for this file. Paths are 
            assumed to have the format:
              - Metadata: "{root}/{fname}/{fname}_{recv}_vis_data"
              - Data:     "{root}/{fname}/{fname}/{fname}_sdp_l0.full.rdb"
        
        bad_ants : list, optional
            List of bad antennas for this file. This will be a list of strings, 
            e.g. ['m012', 'm014'].
            
        diode : dict, optional
            Noise diode settings. Required fields are:
                - 'offset': Timestamp index when diode starts to fire. ('nd_1a0')
                - 'jump_limit': ('f')
                - 'time_range': slice into t_line
        """
        # Validate filename
        filename = valid_filename(filename)
        
        # Check diode kwarg is valid and fill-in default values if not specified
        if diode is None: diode = {}
        if 'offset' not in diode.keys():
            diode['offset'] = None
        if 'jump_limit' not in diode.keys():
            diode['jump_limit'] = 10
        if 'time_range' not in diode.keys():
            diode['time_range'] = (None, None)
        
        # Add entries to info listing
        self.files[filename] = {
            'root':         root,
            'target':       target,
            'diode':        diode,
            'bad_ants':     bad_ants
        }
    
    
