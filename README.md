# JWST_VHS1256b_Reduction
NIRSpec reduction from The JWST Early Release Science Program for Direct Observations of Exoplanetary Systems II: A 1 to 20 Micron Spectrum of the Planetary-Mass Companion VHS 1256-1257 b

The code here was made from examples from Polychronis Patapis and STSci notebooks. 

The code runs individual dithers through the JWST pipeline for observations VHS 1256b and an A-star. The extracted A-star spectra are used to rescale the flux calibration for the VHS 1256 b data.

Caveats:
The accepted version of this spectrum uses jwst pipeline version 1.7.2. Things are still updating in the pipeline so more optimal reductions of the spectrum are very possible. There are lingering issues such as wobbles in the spectrum that may arise from the cube building step when running individual dithers.
