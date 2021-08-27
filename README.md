# SyntheticOTS
OTS scattering is a technique that can be used to determine electron denisty, ion and electron temperatures and velocity in laboratory plasmas. This code generated synthetic OTS spectral density functions for a givens et pf plasma parameters.

Use the function thompsonScatteringCgsV3.

Returns TS Spectral Density Function and the output convolved with the
spectrometer response (gaussian function with s.d. = 2 angstorm)

This code is based on TS formulas derived by Sheffield et al. 2010 

Example Usage:
[S,shift,u,v,w] = thompsonScatteringCgsV3(5e17,6,12,70,70,0,45,532,4)

Author: Datta, Rishabh. 2021. MIT 
