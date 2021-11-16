%% ----- CALCULATE THE REFLECTANCE FUNCTION -----



% By Andrew J. Buggee
%% ------ Read input settings and output data from uv_spec -----

function [R,R_lambda] = reflectanceFunction(inputSettings,ds)

% Geometry values from input Settings -
mu = inputSettings{2}; % cosine of the viewing zenith angle
phi = inputSettings{3}; % sensor aziumuth angle
sza = inputSettings{4}; % solar zenith angle
phi0 = inputSettings{5}; % solar azimuth angle
sensorAlt = inputSettings{6}; % sensor altitude in km
source = inputSettings{7}; % - mW/(m^2 nm) - source irradiance

% radiative transfer solutions
wavelength = ds.wavelength;
irrad0 = source(:,2); % - mW/(m^2 nm) -  source irradiance


% a few other constants calculated from the inputs
mu0 = cosd(sza); % cosine of the solar zenith angle
geomSets = length(mu)*length(phi);

% introduce the spectral responces function
specRep = ones(size(wavelength)); % perfect spectral response

%% ----- CALCULATE REFLECTANCE FUNCTION -----

% reflectance varies with lambda, tau_cloud, r_e, viewing angle, solar
% zenith angle, and the viewing azimuth angle. uvspec can only vary the
% viewing azimuth and viewing zenith angle for a given file. So we will
% create a 3D array that varies these three parameters per file

R_lambda = zeros(length(wavelength),geomSets); % phi (azimuth) changes first, then mu (cos(sza))
R = zeros(length(mu),length(phi));
% if there is a single monochromatic wavelength then we don't need to
% integrate. We simply divide
if length(wavelength)==1
    
    for ii = 1:geomSets
        
        R_lambda(ii) = pi*ds.radiance(ii).rad_umu_phi./(mu0*irrad0); % - 1/sr - reflectance function for monochromatic calculation
        
    end
    
elseif length(wavelength)>1
    
    for ii = 1:geomSets
        
        R_lambda(:,ii) = pi*ds.radiance(ii).rad_umu_phi./(mu0*irrad0); % - 1/sr - reflectance function for monochromatic calculation
        R(ii) = trapz(wavelength,R_lambda(ii).*specRep.*irrad0)./trapz(wavelength,specRep.*irrad0); % - 1/sr - reflectance function over a finite bandwidth
        
    end
    
else
    
    error('Something is wrong with the wavelength vector');
    
end









end