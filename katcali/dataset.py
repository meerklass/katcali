import numpy as np
import json
import pickle
import katdal
from .models import cal_sources
from .utils import valid_filename, file_size


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
            pass
    
    
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
    
    
    def load_metadata(self, filename, verbose=True):
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
        if verbose: print("Metadata file: %3.1f %s" % (fsize, units))
        
        # Load the metadata
        meta = katdal.open(fname)
        self._metadata[filename] = meta
    
    
    def load_data(self, filename, recv, pol, verbose=True):
        """
        Load visibility data from a given file for a given receiver.
        
        Parameters
        ----------
        filename : str
            Name of file to load data for.
        
        recv : str
            Name of receiver to load data for, e.g. 'm006'.
        
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
        recvp = "%s%s" % (recv, pol)
        
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
    
    
    def get_diode_info(self, filename, field=None):
        """
        Get information about the noise diode settings for a particular file.
        
        Parameters
        ----------
        filename : str
            Name of file to find noise diode model for.
        
        field : str, optional
            If specified, return the data for a given setting only. Otherwise, 
            a dictionary of all settings for the file will be returned. 
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
    
    
