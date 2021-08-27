# SyntheticOTS
OTS scattering is a technique that can be used to determine electron denisty, ion and electron temperatures and velocity in laboratory plasmas. This code generated synthetic OTS spectral density functions for a givens et pf plasma parameters.

Use the function thompsonScatteringCgsV3.

Returns TS Spectral Density Function and the output convolved with the
spectrometer response (gaussian function with s.d. = 2 angstorm)
This code is based on TS formulas derived by Sheffield et al. 2010 
%%%%%%%%%
% INPUTS:
% ne = Electron Denisty [cm^-3]
% z = Ionization [-]
% a = Ion mass / proton mass[-]
% te = Electron Temp [eV]
% ti = Ion DenTemp [eV]
% u = plasma velocity along the scattering vector [cm/s]
% th = scattering angle, degrees
% lam_i = incident wavelength, [nm]
% rng = range of wavelength shift, [angstrong]
% OUTPUTS:
% S = spectral density function
% shift = wavelengthshift in range (-rng,rng), Angstrom
% u = spectrometer response
% v = Normalized Spectral Denisty Function S / max(S)
% w = Convolved and Normalized spectral Denisty function
%%%%%%%%%
% Example Usage:
% [S,shift,u,v,w] = thompsonScatteringCgsV3(5e17,6,12,70,70,0,45,532,4)
%%%%%%%%%
% Author: Datta, Rishabh. 2021. MIT 
