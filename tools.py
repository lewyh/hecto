#from scipy.stats import nanmedian
from scipy import optimize
from astropy.io import fits
from pylab import plot, nan
import profiles
import numpy
import glob
import spcl
import os

def fit_continuum(spectrum, wavelengths, order=1):
    """Fits polynomial of user defined order to the supplied spectrum.
    
    Required arguments:
            spectrum        : One-dimensional array containing spectrum
            wavelengths     : One-dimensional array containing wavelengths 
                              corresponding to spectrum.
    """
    fit = numpy.polyfit(wavelengths, spectrum, order)
    poly = numpy.poly1d(fit)
    return poly

def lineheight(spectrum, line_id='OH8399'):
    """docstring for lineheight"""
    # Define linedict, which contains line identifiers as keys, and 
    # [minimum wavelength, maximum wavelength, central_wavelength, line width]
    # as values (minimum and maximum wavelengths define the line and a 
    # sensible amount of continuum on either side.
    linedict = { 'OH8399' : [8390.0,8410.0,8399.17,10.0] }
    line = linedict[line_id]

    spectrum.mask_wl(numpy.array(([line[0],line[1]])),line_id)

    spectrum.mask_wl(numpy.array(([line[0], line[2]-(line[3]/2)],
                                  [line[2]+(line[3]/2),line[1]])), 
                                  line_id+' continuum')
    spectrum.sigma_clipping(line_id+' continuum', name=line_id+' clipped',
                            sigma_high=2, iterations=10)

    fit = fit_continuum(spectrum.spectra[line_id+' clipped'][0],
                        spectrum.spectra[line_id+' clipped'][1], order=1)

    lineprofile = spectrum.spectra[line_id][0] - \
                  fit(spectrum.spectra[line_id][1])

    opt = optimize.leastsq(profiles.errfunc, [line[2],line[3],1],
                           args=(spectrum.spectra[line_id][1],
                                 lineprofile, profiles.gaussian))

    return fit[0][2]

    pass

def redleak_correct(data, critical_slope=0.8):
    """Function which reads in spectrum, fits slope to continuum above 8600A.
    If slope is greatr than some critical value, then assume that red leak 
    needs to be fixed. Accomplish this by fitting higher order polynomial to
    the continuum in this range, and correct spectrum so continuum is flat, 
    taking correct value to be the median count level between 8556-8577 A.
    
    Requires arguments:
            spectrum        : One-dimensional array containing spectrum.
            wavelengths     : One-dimensional array containing wavelength 
                              range (obviously should be same extent as the
                              spectrum array).
    Optional agruments:
            critical_slope  : Value for slope above which red leak correction 
                              is to be applied.
    """
    wavelengths = data.spectra['original'][1]
    spectrum = data.spectra['original'][0]
    header = data.header
    x = [8510,8530]
    mask = numpy.where((wavelengths>x[0])&(wavelengths<x[1]))
    lims = [[8555,8580],[8720,8755],[8782,8788],[8796,8806],[8814,8822],
            [8854,8860],[8872,8882],[8890,8890],[8924,8938],[8947,8954],
            [8965,8980],[9005,9035],[9040,9045]]
    for x in lims:
        masktemp = numpy.where((wavelengths>x[0])&(wavelengths<x[1]))
        mask = numpy.append(mask, masktemp)

    fit, slope = fit_continuum(spectrum[mask],wavelengths[mask])

    if slope>critical_slope:
        fit, slope = fit_continuum(spectrum[mask],wavelengths[mask], order=3)
        continuum = spectrum[numpy.where((wavelengths>8556)&\
                                         (wavelengths<8577))]
        flat_value = numpy.median(continuum)
        fixmask = numpy.where(wavelengths>8570)[0]
        redval = [0]*len(fixmask)
        for i in range(len(fixmask)):
            redval[i] = fit(wavelengths[fixmask[i]])
        subtract = redval-flat_value
        corrected = numpy.array(spectrum[fixmask]-subtract,
                                dtype=numpy.float32)
        spectrum[fixmask]=corrected
        for i in range(len(fit.coeffs)):
            coeffi = str(fit.coeffs[i])
            header.update('REDLEAK'+str(i+1),coeffi,
                          comment="Coefficient of order "+\
                                  str(i)+" of polynomial fit")
        header.update('REDLEAKV',flat_value,comment="Value for flat cont.")
    
    return spectrum, header
