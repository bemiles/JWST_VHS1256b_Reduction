from astropy.table import Table
from specutils import Spectrum1D
from header_function import *
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np

fits_data_dir = '/Users/jolie/PycharmProjects/JWST_VHS1256b_Reduction/write_spectra/'  #Where you want to write fits files to

#Header information about both NIRSpec and MIRI Spectra
JWST_header_dict = {
    'RA': '194.007637',
    'Dec': '-12.957692',
    'generated_history': 'Data Reduced by JWST pipeline version 1.7.2, CRDS version 11.16.12, CRDS context jwst_0977.pmap',
    'object_name' : 'VHS1256b',
    'telescope': 'JWST',
    'comment': None
     }


#################Generate NIRSpec Files ##################

# Input your NIRSpec wavelength, flux, and error data here
# Be sure to include units
NIRspec_wavelength = np.arange(0, 5, 0.5) * u.um
NIRspec_flux = [0.15, .15, .07,.05,.27,.022, .06, .037, .02, .2 ] * u.Jy
NIRspec_error = [0.0014, .006, .0088, .009, .0097, .0098, .007, .0013, .0012, .009] * u.Jy

NIR_table = Table({'wavelength': NIRspec_wavelength,
                   'flux': NIRspec_flux,
                   'flux_uncertainty': NIRspec_error})

NIRSpec_filename = 'VHS1256b_NIRSpec'
fits_filename_NIR = fits_data_dir + NIRSpec_filename + '.fits'

# Header Info specific to NIRSpec spectra
NIRSpec_header_dict = {
    'instrument': 'NIRSpec',
    'bandpass' : 'NIR'
     }


NIRSpec_dict = {**JWST_header_dict, **NIRSpec_header_dict} # combine dictionaries to create headers

# Turns dictionary ino properly formatted FITS header
header = compile_header(NIR_table['wavelength'], **NIRSpec_dict)

# Make a dedicated HDU for the header and spectrum
hdu0 = fits.PrimaryHDU(header=header)
hdu1 = fits.BinTableHDU(data=NIR_table)
hdu1.header['EXTNAME'] = 'SPECTRUM'
hdu1.header.set('OBJECT', NIRSpec_dict['object_name'], 'Object Name')

# Combine the two HDUs into a multi-extension FITS (MEF)
spectrum_mef = fits.HDUList([hdu0, hdu1])  # hdu0 is header and hdu1 is data

# Write out the MEF to a file
spectrum_mef.writeto(fits_filename_NIR, output_verify="exception")

# This is extra code which demonstrates how to read the file back in using specutils
spec1d_NIR = Spectrum1D.read('/Users/jolie/PycharmProjects/JWST_VHS1256b_Reduction/write_spectra/VHS1256b_NIRSpec.fits', format='tabular-fits')
header_NIR = fits.getheader(fits_filename_NIR)
name_NIR = header_NIR['OBJECT'] # get object name from the header

# Plot spectra
ax = plt.subplots()[1]
ax.plot(spec1d_NIR.spectral_axis, spec1d_NIR.flux)
ax.errorbar(spec1d_NIR.spectral_axis.value, spec1d_NIR.flux.value, yerr=spec1d_NIR.uncertainty.array, fmt='-')
ax.set_xlabel(f"Dispersion ({spec1d_NIR.spectral_axis.unit})")
ax.set_ylabel(f"Flux ({spec1d_NIR.flux.unit})")
plt.title(f"{name_NIR} {header_NIR['TELESCOP']} {header_NIR['INSTRUME']}") # pulls spectra info from the FITS header
plt.savefig(fits_data_dir + NIRSpec_filename + '.png')


#######Generate MIRI Files #######################

# Input your NIRSpec wavelength, flux, and error data here
# Be sure to include units
MIRI_wavelength = np.arange(5, 15, 1) * u.um # and 5-20 microns (MIRI)
MIRI_flux = [0.15, .15, .07,.05,.27,.022, .06, .037, .02, .2 ]* u.Jy
MIRI_error = [0.0014, .006, .0088, .009, .0097, .0098, .007, .0013, .0012, .009] * u.Jy


MIRI_table = Table({'wavelength': MIRI_wavelength,
                           'flux': MIRI_flux,
                           'flux_uncertainty': MIRI_error})

MIRI_filename = 'VHS1256b_MIRI'
fits_filename_MIRI = fits_data_dir + MIRI_filename + '.fits'

# Header Info specific to MIRI spectra
MIRI_header_dict = {
    'instrument': 'MIRI',
    'bandpass' : 'MIR'
     }

MIRI_dict = {**JWST_header_dict, **MIRI_header_dict}  # combine dictionaries to create headers

# Turns dictionary ino properly formatted FITS header
header = compile_header(MIRI_table['wavelength'], **MIRI_dict)

# Make a dedicated HDU for the header and spectrum
hdu0 = fits.PrimaryHDU(header=header)
hdu1 = fits.BinTableHDU(data=MIRI_table)
hdu1.header['EXTNAME'] = 'SPECTRUM'
hdu1.header.set('OBJECT', MIRI_dict ['object_name'], 'Object Name')

# Combine the two HDUs into a multi-extension FITS (MEF)
spectrum_mef = fits.HDUList([hdu0, hdu1])  # hdu0 is header and hdu1 is data

# Write out the MEF to a file
spectrum_mef.writeto(fits_filename_MIRI, output_verify="exception")


# This is extra code which demonstrates how to read the file back in using specutils
spec1d_MIRI = Spectrum1D.read('/Users/jolie/PycharmProjects/JWST_VHS1256b_Reduction/write_spectra/VHS1256b_MIRI.fits', format='tabular-fits')
header_MIRI = fits.getheader(fits_filename_MIRI)
name_MIRI = header_MIRI['OBJECT'] # get object name from header

# Plot spectra
ax = plt.subplots()[1]
ax.plot(spec1d_MIRI.spectral_axis, spec1d_MIRI.flux)
ax.errorbar(spec1d_MIRI.spectral_axis.value, spec1d_MIRI.flux.value, yerr=spec1d_MIRI.uncertainty.array, fmt='-')
ax.set_xlabel(f"Dispersion ({spec1d_MIRI.spectral_axis.unit})")
ax.set_ylabel(f"Flux ({spec1d_MIRI.flux.unit})")
plt.title(f"{name_MIRI} {header_MIRI['TELESCOP']} {header_MIRI['INSTRUME']}") #pulls spectra info from the FITS header
plt.savefig(fits_data_dir + MIRI_filename + '.png')


