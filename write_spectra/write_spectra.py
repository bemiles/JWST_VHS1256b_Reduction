import os
from urllib.parse import unquote
import astropy.io.fits as fits
from astropy.table import Table
#from specutils import Spectrum1D
from header_function import *
import astropy.units as u


fits_data_dir = '/Users/jolie/gitlocation/SIMPLE-db/scripts/test_scripts/'  # set outside the dictionary


#first generate wavelength, flux, and error arrays
#FAKE DATA
NIRspec_wavelength = [.6,.7,.8,.9,1,1.1,1.2,1.4,1.4,1.5]  #1-5 microns (NIRSpec)
NIRspec_flux = [-0.15, .15, .07,.05,.27,.022, .06, .037, .02, .2 ]
NIRspec_error = [0.14, .06, .088, .09, .097, .098, .07, .013, .012, .09]

NIR_table = create_NIR_table(NIRspec_wavelength, NIRspec_flux,NIRspec_error ) #table of data


MIRI_wavelength = [5.6,5.7,5.8,5.9,5,6.1,6.2,6.4,6.4,6.5]   # and 5-20 microns (MIRI)
MIRI_flux = [-0.15, .15, .07,.05,.27,.022, .06, .037, .02, .2 ]
MIRI_error = [0.14, .06, .088, .09, .097, .098, .07, .013, .012, .09]

MIRI_table = create_MIRI_table(MIRI_wavelength, MIRI_flux,MIRI_error ) #table of data

JWST_header_dict = {
    # Information about new spectra
    'RA': '',
    'Dec': '',
    'generated_history': 'Data Reduced by JWST pipeline version 1.7.2, CRDS version 11.16.12, CRDS context jwst_0977.pmap',
    'object_name' : 'VHS1256b',
    'telescope': 'JWST',
    'comment': None
     }


NIRSpec_header_dict = {
    'instrument': 'NIRSpec',
    'bandpass' : 'NIR',
    'generated_filename': 'VHS1256b_NIRSpec'
     }

NIRSpec_dict = {**JWST_header_dict, **NIRSpec_header_dict}



MIRI_header_dict = {
    'instrument': 'MIRI',
    'bandpass' : 'MIR',
    'generated_filename': 'VHS1256b_MIRI'
     }

MIRI_dict = {**JWST_header_dict, **MIRI_header_dict}



#######GENERATE NIRSpec Fits####################################
object_name = unquote(NIRSpec_dict['object_name'])
NIRSpec_dict ['object_name'] = object_name

header = compile_header(NIR_table['wavelength'], **NIRSpec_dict)

spectrum_data_out = Table({'wavelength': NIR_table['wavelength'], 'flux': NIR_table['flux'], 'flux_uncertainty': NIR_table['error'] })


# Make the HDUs
hdu1 = fits.BinTableHDU(data=spectrum_data_out)
hdu1.header['EXTNAME'] = 'SPECTRUM'
hdu1.header.set('OBJECT', object_name, 'Object Name')
hdu0 = fits.PrimaryHDU(header=header)

# Write the MEF with the header and the data
spectrum_mef = fits.HDUList([hdu0, hdu1])  # hdu0 is header and hdu1 is data

fits_filename = fits_data_dir + NIRSpec_dict['generated_filename'] + '.fits'
try:
    spectrum_mef.writeto(fits_filename, overwrite=False, output_verify="exception")
except:
    raise



#######GENERATE MIRI Fits####################################
object_name = unquote(MIRI_dict['object_name'])
MIRI_dict ['object_name'] = object_name

header = compile_header(MIRI_table['wavelength'], **MIRI_dict)

spectrum_data_out = Table({'wavelength': MIRI_table['wavelength'], 'flux': MIRI_table['flux'], 'flux_uncertainty': MIRI_table['error'] })


# Make the HDUs
hdu1 = fits.BinTableHDU(data=spectrum_data_out)
hdu1.header['EXTNAME'] = 'SPECTRUM'
hdu1.header.set('OBJECT', object_name, 'Object Name')
hdu0 = fits.PrimaryHDU(header=header)

# Write the MEF with the header and the data
spectrum_mef = fits.HDUList([hdu0, hdu1])  # hdu0 is header and hdu1 is data

fits_filename = fits_data_dir + MIRI_dict['generated_filename'] + '.fits'
try:
    spectrum_mef.writeto(fits_filename, overwrite=False, output_verify="exception")
except:
    raise






"""spec1d = Spectrum1D.read(fits_filename, format='tabular-fits')
header = fits.getheader(fits_filename)
name = header['OBJECT']


ax = plt.subplots()[1]
# ax.plot(spec1d.spectral_axis, spec1d.flux)
ax.errorbar(spec1d.spectral_axis.value, spec1d.flux.value, yerr=spec1d.uncertainty.array, fmt='-')
ax.set_xlabel(f"Dispersion ({spec1d.spectral_axis.unit})")
ax.set_ylabel(f"Flux ({spec1d.flux.unit})")
plt.title(f"{name} {header['TELESCOP']} {header['INSTRUME']}")
plt.savefig(NIRSpec_dict['fits_data_dir'] + NIRSpec_dict['generated_filename'] + '.png')

### ^ rn these are outside of dictionary

"""

