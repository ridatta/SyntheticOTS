function [S,shift,u,v,w] = thompsonScatteringCgsV3(ne,z,a,te,ti,u,th,lam_i,rng)
% Returns TS Spectral Density Function and the output convolved with the
% spectrometer response (gaussian function with s.d. = 2 angstorm)
% This code is based on TS formulas derived by Sheffield et al. 2010 
%%%%%%%%%
% INPUTS:
% ne = Electron Denisty [cm^-3]
% z = Ionization [-]
% a = Ion mass / proton mass[-]
% te = Electron Temp [eV]
% ti = Ion DenTemp [eV]
% u = plasma velocity along teh scattering vector [cm/s]
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
%%%%%%%%%

global me mi Te Ti Z lamD kb A V 

% Load constants
load physicalConstants-cgs.mat c me kb mp 

% Parameters
n0 = ne; % [cm^-3], unperturbed electron denisty
Z = z; % ionization 
A = a; % ion atomic mass
Te = te; % Electron Temp, eV
Ti = ti; % Ion Temp, eV
rng = rng * 1; % [Arngstorm]
V = u; % Velocity along k, (k dot V / k) [cm/s]
% Geometrical parameters
theta = deg2rad(th); % scattering angle, rad
lambda_i = lam_i * 1e-9 * 100; % incident wavelength, [cm]

% calculated parameters
mi = A * mp; % ion mass, [g]
ki = 2 * pi / lambda_i; % incident waveumber [cm^-1]
FP = FundamentalPlasma(); % Load the fundamental plasma object
w_pe = FP.getElecPlasmaFreq(n0); % rad/s
lamD = sqrt(2) * FP.getDebyeLen(n0,Te); % [cm^-1]

shift = linspace(-rng,rng,100*7.5) * 1e-10 * 100; % wavelength shift, [cm]
S = zeros(numel(shift),1); 
% default number of points is 750

for ii = 1:length(shift)
delL = shift(ii); % wavelength shift, [cm]
lambda_s = lambda_i + delL; % scattered wavelength, [cm]
% wave vectors
ks = 2 * pi / lambda_s; % scttered wavenumber [cm^-1]
k = sqrt(ki^2 + ks^2 - 2 * ks * ki * cos(theta)); % scattering vector, [cm^-1]
% frequency
wi = sqrt(c^2 * ki^2 + w_pe^2); % incident freq, [s^-1]
ws = sqrt(c^2 * ks^2 + w_pe^2) + V * k; % scattered freq, [s^-1]
w = ws - wi; % scattering freq, [s^-1]

S(ii) = getSpectralDensity(k,w); % form factor
end

% Convolve with gaussian function
v = S/max(S); % normalized spectral denisty
u = 0.25 * gaussianFn(shift,0,0.2e-10*100); % spectrometer response, gaussian function 
w = conv(v,u','same'); % convolve

% % Plot spectral denisty
% figure
% plot(shift * 10^8, v, 'k-', 'LineWidth', 3, 'DisplayName','Spectral Density'); hold on;
% xlabel('$\lambda_i - \lambda_s, \AA$','Interpreter','latex'); 
% ylabel('$\hat{S}({\bf k},\omega)$','Interpreter','latex'); 
% formatPlots(); 
% legend('location','northwest');
% plot(shift * 10^8, u, 'k--','LineWidth', 3,'handleVisibility','off'); hold on;
% plot(shift * 10^8, w/max(w), 'r-','LineWidth', 3, 'DisplayName','Convolved'); hold on;
% fname = ['n0', num2str(n0), '-', ...
%             'Z', num2str(Z), '-',...
%              'A', num2str(A), '-',...
%              'Te', num2str(Te), '-',...
%              'Ti', num2str(Ti), '-',...
%              'theta', num2str(theta), '-', ...
%              'wl',num2str(lambda_i*1e9/100)];     
% title(['$\lambda_i = ' num2str(lambda_i/100 *1e9) ' \,  nm, ' ...
%     '\theta = ' num2str(rad2deg(theta)) '^o $'], 'Interpreter','latex'); 
% save([fname, '.mat']);          
% saveas(gcf, [fname '.png']); 
end

function out = gaussianFn(x,mu,sig)
% Returns the normalized gaussian fn
% mu = mean 
% sig = standard deviation
y = 1 / (sig * sqrt(2 * pi)) * exp(-0.5 * ((x - mu)/sig).^2); 
out = y / max(y); 
end

function out = getSpectralDensity(k,w) % Caluclates the output at a given wavelength
global Te Ti Z A
addpath(checkDir('/Users/rishabhdatta/Dropbox (MIT)/PUFFIN/Codes/Other/PlasmaFormulary\'));
FP = FundamentalPlasma(); 
vthe = sqrt(2) * FP.getElecThermalVel(Te); 
vthi = sqrt(2) * FP.getIonThermalVel(Ti,A); 

chi_e = getChi(k,w,vthe,1); % electron suceptibility
chi_i = getChi(k,w,vthi,Z); % ion suceptibility
ep = chi_e + chi_i + 1; % dielectric fn

fe0 = @(v) FP.getMaxwellDist1D(v,vthe); % electron distribution
fi0 = @(v) FP.getMaxwellDist1D(v,vthi); % ion distribution

out = 2 * pi / k *  (abs(1 - chi_e / ep))^2 * fe0(w/k) + ...
    2 * pi * Z / k * (abs(chi_e / ep))^2 * fi0(w/k); 
end

function out = getChi(k,w,vth,zz) % Get succeptibility
global lamD Te Ti 
W = @(z0) 1 + 1i * sqrt(pi) * z0 * exp(-z0^2) - 2 * z0 * dawson(z0); 
out = 1 / (lamD * k)^2 * (zz * Te / Ti) * W(w/(k * vth)); 
end


