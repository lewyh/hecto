from astropy.io import fits
import numpy as np
from tools import *
import math
import os
from pylab import *
numpy.set_printoptions(threshold=3)


def INDEX2WAVE(INDEX, CRVAL1, CD1):
    WAVE = CRVAL1+(INDEX*CD1)
    return WAVE

def WAVE2INDEX(WAVE, CRVAL1, CD1):
    INDEX = (WAVE-CRVAL1)/CD1
    return INDEX 

class Spectrum:
    def __init__(self, filename, ext=0, start=6532.0, end=9045):
        """Reads a spectrum into an array, can deal with multi-extension
           Hectospec files.

        Required arguments:
                filename    : Path to file containing spectrum of interest

        Optional arguments:
                spectype    : Keyword specifying what kind of spectrum is
                              contained in the file - usually 'object' or 
                              'sky'. Defaults to 'object'.
                ext         : Which extension should be opened. For 
                              Hectospec files, this is usually 0 (default).
                ap          : Which aperture inside the extension contains
                              the spectrum of interest. For Hectospec reduced
                              files,
                                    0 = object - (scaled sky)
                                    1 = raw spectrum
                                    2 = averaged sky
                start/end   : Wavelength values defining range of interest.
                              The entire spectrum will be loaded, but 
                              mask['true'] will define this region.

        Creates:
            self.wavelengths        : array of wavelengths in start/end range
            self.all_wavelengths    : array of all wavelengths in file
            self.spectrum           : array of spectrum in start/end range
            self.all_spectrum       : array of entire spectrum in file
            self.CRVAL1             : CRVAL1,
            self.CD1                : CD1
            self.mask['true']       : mask of shape spectrum, all == True.
            self.spectra['original']: dict which makes it easier to saved
                                      multiple versions of the spectra 
                                      (e.g. normalized etc.).
        """
        self.header = fits.getheader(filename, 0)
        if '/' in filename:
            filen = filename.split('/')[-1]
            self.directory = filename.rstrip(filen)
        else:
            self.directory = './'
        CRVAL1, CD1 = self.header['CRVAL1'], self.header['CD1_1']

        self.raw_file = fits.open(filename)

        if self.header['INSTRUME']=='hectospec':
            self.units = self.header['WAT1_001'].split()[-1].split('=')[-1]
            self.target = self.header['POINTING']+'/'+\
                          self.header['APERTURE']+'.'+self.header['CATOBJ']+\
                          ".ms.fits"
            if self.header['CATOBJ']=='sky':
                self.spectype = 'sky'
            else:
                self.spectype = 'object'
            # Work only with Hectospec data for now. Read in the first fits 
            # extension (the only one used by Hectospec). If the length of the 
            # array in this extension is greater than one, then import the 
            # array number specified by the user. Otherwise open the only 
            # existing array.
            self.data = fits.getdata(filename, 0)
            if len(self.data)>1:
                onedspectrum = self.data[ext][0]
            else:
                onedspectrum = self.data[0][0]
        # Calculate the maximum wavelength.
        WAVELMAX = INDEX2WAVE(len(onedspectrum), CRVAL1, CD1)
        # Calculate array of wavelengths to simplify plotting
        wavelengths = np.arange(CRVAL1, WAVELMAX, CD1)
            
        mask = (wavelengths > start)&(wavelengths < end)
        start = WAVE2INDEX(start, CRVAL1, CD1)
        end = WAVE2INDEX(end, CRVAL1, CD1)

        self.wavelengths = wavelengths
        self.spectrum = onedspectrum
        self.CRVAL1 = CRVAL1
        self.CD1 = CD1
        self.masks = {'true':mask}
        self.spectra = {'original':[onedspectrum, wavelengths]}
        self.filename = filename
        self.ext = ext

    def add_spectrum(self, spectrum, wavelengths, spec_id):
        """docstring for add_spectrum"""
        
       # self.spectra['spec_id']=
        pass

    def _use_spectrum(self, spectrum):
        try:
            spectrum = self.spectra[spectrum]
        except:
            print "You have tried to use a spectrum that does't exist."
            print "Check that you are referencing a valid spectrum."
            print "You tried to use spectrum key \'"+spectrum+"\'"
            print "The spectrum keys available are: "+str(self.spectra.keys())
            return 0
        return spectrum[0], spectrum[1]

    def _use_mask(self, mask):
        try:
            mask = self.masks[mask]
        except:
            print "You have tried to use a mask that does't exist."
            print "Check that you have created the mask."
            print "You tried to use mask \'"+mask+"\'"
            print "The masks available are: "+str(self.masks.keys())
            return 0
        return mask

    def plot(self, spec_id='original', mask='true', 
             fig=None, axes=None, subplot=111, v_offset=0, titles=True):
        spectrum, wavelengths = self._use_spectrum(spec_id)
        if (type(spectrum)==int):
            return
        else:
            if get_fignums()==[] or fig=='new':
                fig = figure()
                axes = fig.add_subplot(subplot)
            xlim(min(wavelengths),
                 max(wavelengths))
            xlabel('Wavelength ('+self.units+')')
            ylabel('Counts')
            if titles==True:
                suptitle(self.target,size=14)
                title("("+self.spectype+" "+self.header['RA']+" "+\
                      self.header['DEC']+" EXPTIME:"+\
                      str(self.header['EXPTIME'])+"s)",size=12)
            plot(wavelengths,spectrum+v_offset, 
                 label=self.target)
        return axes

    def mask_wl(self, wavelength_ranges, name, 
                spec_id='original', invert=False):
        """Accepts spectrum from read(). Returns mask based on provided
           wavelength ranges.

        Requires arguments:
                range_name        : unique identifier so that mask can
                                    be saved and referenced more easily
                wavelength_ranges : either an array of length two if only one
                                    range is needed, or a numpy array of
                                    shape (N,2) if multiple ranges are desired
        Optional arguments:
                invert            : returns inverted mask w.r.t. range(s) 
                (default = False)   supplied - useful for masks that simply 
                                    exclude a region (e.g. remove a line). 
        """
        spectrum, wavelengths = self._use_spectrum(spec_id)
        mask = wavelengths>0
        # Checks shape of supplied wavelength range array in order to
        # determine whether one or more ranges have been supplied. Can't
        # loop over a single range so treat this case differently.
        if wavelength_ranges.shape == (2,):
            temp_mask = (wavelengths>wavelength_ranges[0])&\
                        (wavelengths<wavelength_ranges[1])
            final_mask = temp_mask&mask
        else:
            # Create flag to differentiate between first item in array
            # and the rest - in order to initialize the final_mask.
            flag = 0
            for each_pair in wavelength_ranges:
                # For each wavelength range, create a mask for that range
                temp_mask = (wavelengths>each_pair[0])&\
                            (wavelengths<each_pair[1])
                # Combine it with the mask supplied to wl_mask through
                # spectrum[4] in case any values have already been masked
                temp_mask2 = temp_mask&mask
                if flag==0:
                    # If this is the first wavelength range, create
                    # final_mask, which will be combined with the subsequent
                    # ranges' masks.
                    final_mask = temp_mask2
                    flag=1
                else:
                    # Combine the wavelength ranges' masks with an OR
                    # operator, marking all values in any of the desired
                    # ranges as True.
                    final_mask = temp_mask2|final_mask
        # If an inverted mask has been requested, then inverse final_mask.
        if invert==True:
            final_mask=~final_mask

        # Saves the mask as 'name' in the self.masks dict. Also saves the 
        # masked [spectrum, wavelength] as a list in the self.spectra dict.
        self.masks[name]=final_mask
        self.spectra[name]=[spectrum[final_mask], wavelengths[final_mask]]
        return

    def redleak(self, spec_id='original'):
        redleak_returned = redleak_correct(self)
        self.spectra['corrected']=redleak_returned[0]
        self.header = redleak_returned[1]
        return

    def sigma_clipping(self, spec_id='original', name='clipped',
                       sigma_high=3, iterations=3):
        spectrum, wavelengths = self._use_spectrum(spec_id)
        i=0
        while i<iterations:
            median = numpy.median(spectrum)
            stdev = numpy.std(spectrum)
            diff = abs(spectrum-median)
            if i==0:
                mask = diff<(sigma_high*stdev)
            else:
                mask = mask & (diff<(sigma_high*stdev))
            i+=1
        self.masks[name] = mask
        self.spectra[name] = [spectrum[mask], wavelengths[mask]]
        return mask

