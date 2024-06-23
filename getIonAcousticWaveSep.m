function out = getIonAcousticWaveSep(n0,Te,Ti,Z,A,th,lam_i)
% n0 = elecrron denisty [cm^-3]
% Te, Ti = [eV]
% Z = avg. ionization [-]
% A = atomic mass [-]
% th = scattering ang. [deg]
%  lam_i = wavelength of light = [nm]

addpath('/Users/rishabhdatta/Dropbox (MIT)/PUFFIN/Codes/Other/PlasmaFormulary')

F = FundamentalPlasma();

lam_d = F.getDebyeLen(n0,Te); % cm


load physicalConstants-SI.mat mp E0

lam_i = lam_i * 1e-9; % m
c = 3e8; % m/s
th = deg2rad(th); % rad
Te = E0 * Te; % J
Ti = E0 * Ti; % J
mi = A * mp; % kg
lam_d = cgs2SI(lam_d,'Length'); % m

k = 2 * pi / lam_i; % m^-1

out = lam_i * 4/c * sin(th/2) *...
    sqrt( Te/mi * ( Z / (1 + k^2*lam_d^2) + 3 * Ti / Te));

out = out * 1e10; % Anstrom

end