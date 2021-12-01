%% ----- READ IN MODIS DATA -----

% this function will read in L1B MODIS data as .hdf files
% will produce all necessary information into a cell array, and the data
% set into a structure

% there are many different fields of data one could read from a MODIS .hdf


% ---- Description of the different fields within the HDF File -----

%   (1) radiance - units of W/m^2/micron/sr
%   (2) reflectance - units of 1/sr



% By Andrew J. Buggee
%%

function [EV] = readMODIS_L1B_data(fileName)
%% ---- Read in Conversion Scales and Offsets -----

% retreive the scales and offsets for radiance and reflectance. Each will
% be a vector with the number of entities equal to the number of bands.

radianceScales = hdfread(fileName,'radiance_scales'); %  output will be a cell array
radianceScales = radianceScales{1};

radianceOffsets = hdfread(fileName,'radiance_offsets'); % output will be a cell array
radianceOffsets = radianceOffsets{1};

reflectanceScales = hdfread(fileName,'reflectance_scales');
reflectanceScales = reflectanceScales{1};

reflectanceOffsets = hdfread(fileName,'reflectance_offsets');
reflectanceOffsets = reflectanceOffsets{1};

%% --- Read in a data set at a specifc resolution ---

% Check to see what product we are reading in

if strcmp(fileName(6:8),'QKM')
    % then we are reading in calibrated Earth View data at 250m resolution
    
    
    uncertainty = 'EV_250_RefSB_Uncert_Indexes';
    
    earthView250 = hdfread(fileName,'EV_250_RefSB');
    
    m250_scaledIntegers = bandAcrossAlong2AcrossAlongBand(earthView250);
    
    EV.m250.radiance = scalesOffsets2Matrix(m250_scaledIntegers,radianceScales,radianceOffsets);
    EV.m250.reflectance = scalesOffsets2Matrix(m250_scaledIntegers,reflectanceScales,reflectanceOffsets);
    
    % the 250 meter resolution data covers the first two MODIS bands
    EV.m250.bands = readMODISbands([1,2]);
    
    
elseif strcmp(fileName(6:8),'HKM')
    % then we are reading in calibrated Earth View data at 500m resolution,
    % including the 250m resolution bands aggregated to 500m resolution
    
    
    uncertainty = 'EV_500_RefSB_Uncert_Indexes';
    
    earthView250 = hdfread(fileName,'EV_250_Aggr500_RefSB');
    earthView500 = hdfread(fileName,'EV_500_RefSB');
    
    m250_scaledIntegers = bandAcrossAlong2AcrossAlongBand(earthView250);
    m500_scaledIntegers = bandAcrossAlong2AcrossAlongBand(earthView500);
    
    
    EV.m250.radiance = scalesOffsets2Matrix(m250_scaledIntegers,radianceScales,radianceOffsets);
    EV.m250.reflectance = scalesOffsets2Matrix(m250_scaledIntegers,reflectanceScales,reflectanceOffsets);
    
    % the 250 meter resolution data covers the first two MODIS bands
    EV.m250.bands = readMODISbands([1,2]);
    
    EV.m500.radiance = scalesOffsets2Matrix(m500_scaledIntegers,radianceScales,radianceOffsets);
    EV.m500.reflectance = scalesOffsets2Matrix(m500_scaledIntegers,reflectanceScales,reflectanceOffsets);
    
    % the 500 meter resolution data covers bands 3 through 7
    EV.m500.bands = readMODISbands([3:7]);
    
    
    
elseif strcmp(fileName(6:8),'1KM')
    % then we are reading in calibrated Earth View data at 1km resolution
    % including the 250m and 500m data resolution data aggregated to 1km
    % resolution
    
    uncertainty = 'EV_1KM_RefSB_Uncert_Indexes';
    
    earthView250 = hdfread(fileName,'EV_250_Aggr1km_RefSB'); % first two modis bands aggregated to 1km resolution
    earthView500 = hdfread(fileName,'EV_500_Aggr1km_RefSB'); % modis bands 3-7 aggregated to 1km resolution
    earthView1000 = hdfread(fileName,'EV_1KM_RefSB'); % modis bands 8 - 36 at 1km resolution
    
    m250_scaledIntegers = bandAcrossAlong2AcrossAlongBand(earthView250);
    m500_scaledIntegers = bandAcrossAlong2AcrossAlongBand(earthView500);
    m1000_scaledIntegers = bandAcrossAlong2AcrossAlongBand(earthView1000);
    
    EV.m250.radiance = scalesOffsets2Matrix(m250_scaledIntegers,radianceScales,radianceOffsets);
    EV.m250.reflectance = scalesOffsets2Matrix(m250_scaledIntegers,reflectanceScales,reflectanceOffsets);
    
    % the 250 meter resolution data covers the first two MODIS bands
    EV.m250.bands = readMODISbands([1,2]);
    
    EV.m500.radiance = scalesOffsets2Matrix(m500_scaledIntegers,radianceScales,radianceOffsets);
    EV.m500.reflectance = scalesOffsets2Matrix(m500_scaledIntegers,reflectanceScales,reflectanceOffsets);
    
    % the 500 meter resolution data covers bands 3 through 7
    EV.m500.bands = readMODISbands([3:7]);
    
    % --- DONT NEED BANDS 8-36 FOR NOW ---
    
%     % retrieve the 1km bands, bands 8-36
%     EV.km1.radiance = scalesOffsets2Matrix(m1000_scaledIntegers,radianceScales,radianceOffsets);
%     EV.km1.reflectance = scalesOffsets2Matrix(m1000_scaledIntegers,reflectanceScales,reflectanceOffsets);
%     
%     % the 1 kilometer resolution reflective bands span 8-16 and 26
%     EV.km1.bands = readMODISbands([8:16,26]);
    
    
elseif strcmp(fileName(6:8),'OBC')
    % then we are reading in the On-Board Calibrator and Engineering Data
    
    
    
    
    
else
    
    error('Filename is not valid! Check to see it is a L1B MODIS data file')
    
end







end





