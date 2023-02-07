import os
from urllib.parse import unquote
import astropy.io.fits as fits
from astropy.table import Table
from matplotlib import pyplot as plt
from specutils import Spectrum1D

#first generate wavelength, flux, and error arrays

NIRspec_wavelength = []  #1-5 microns (NIRSpec)
NIRspec_flux = []
NIRspec_error = []

MIRI_wavelength = []   # and 5-20 microns (MIRI)
MIRI_flux = []
MIRI_error = []


fits_data_dir = '/Users/fits_output_path',  # set outside the dictionary
#'generated_filename': ' '

"""JWST_header_dict = {
    # Information about new spectra
    'RA': '',
    'Dec': '',
    'generated_history': 'Data Reduced by JWST pipeline version 1.7.2, CRDS version 11.16.12, CRDS context jwst_0977.pmap',
    'object_name' : 'VHS1256b',
    'telescope': 'JWST'
     }"""

#NIRSpec_generated_filename: ' '
"""NIRSpec_header_dict = {
    'instrument': 'NIRSpec',
    'bandpass' = 'NIR'
     }"""
#NIRSpec_dict = {**JWST_header_dict, **NIRSpec_header_dict}



#MIRI_generated_filename: ' '
"""MIRI_header_dict = {
    'instrument': 'MIRI',
    'bandpass' = "MIR" 
     }"""

#NMIRI_dict = {**JWST_header_dict, **MIRI_header_dict}




dictionary =

object_name = unquote(spectrum_info_all['object_name'])
spectrum_info_all['object_name'] = object_name

spectrum_path = spectrum_info_all['file_path']
file = os.path.basename(spectrum_path)




header = compile_header(wavelength, **spectrum_info_all)

spectrum_data_out = Table({'wavelength': wavelength, 'flux': flux, 'flux_uncertainty': flux_unc})


# Make the HDUs
hdu1 = fits.BinTableHDU(data=spectrum_data_out)
hdu1.header['EXTNAME'] = 'SPECTRUM'
hdu1.header.set('OBJECT', object_name, 'Object Name')
hdu0 = fits.PrimaryHDU(header=header)

# Write the MEF with the header and the data
spectrum_mef = fits.HDUList([hdu0, hdu1])  # hdu0 is header and hdu1 is data

fits_filename = spectrum_info_all['fits_data_dir'] + spectrum_info_all['generated_filename'] + '.fits'
try:
    spectrum_mef.writeto(fits_filename, overwrite=False, output_verify="exception")
    # TODO: think about overwrite
    logger.info(f'Wrote {fits_filename}')
except:
    raise

spec1d = Spectrum1D.read(fits_filename, format='tabular-fits')
header = fits.getheader(fits_filename)
name = header['OBJECT']


ax = plt.subplots()[1]
# ax.plot(spec1d.spectral_axis, spec1d.flux)
ax.errorbar(spec1d.spectral_axis.value, spec1d.flux.value, yerr=spec1d.uncertainty.array, fmt='-')
ax.set_xlabel(f"Dispersion ({spec1d.spectral_axis.unit})")
ax.set_ylabel(f"Flux ({spec1d.flux.unit})")
plt.title(f"{name} {header['TELESCOP']} {header['INSTRUME']}")
plt.savefig(spectrum_info_all['fits_data_dir'] + spectrum_info_all['generated_filename'] + '.png')





