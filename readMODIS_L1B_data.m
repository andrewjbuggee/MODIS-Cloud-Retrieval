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

% retrieve the info hdf info structure
info = hdfinfo(fileName);



%% --- Read in a data set at a specifc resolution ---

% Check to see what product we are reading in

if strcmp(fileName(6:8),'QKM')
    % then we are reading in calibrated Earth View data at 250m resolution
    
    % --- Radiance scales for the first two bands ---
    radianceScales_250m = info.Vgroup.Vgroup(2).SDS(1).Attributes(6).Value; % output should be a vector with 2 entires for the first two bands
    radianceOffset_250m = info.Vgroup.Vgroup(2).SDS(1).Attributes(7).Value; % output should be a vector with 2 entires for the first two bands
    
    % --- Reflectance scales for the first two bands ---
    reflectanceScales_250m = info.Vgroup.Vgroup(2).SDS(1).Attributes(9).Value; % output should be a vector with 2 entires for the first two bands
    reflectanceOffset_250m = info.Vgroup.Vgroup(2).SDS(1).Attributes(10).Value; % output should be a vector with 2 entires for the first two bands
    
    uncertainty = 'EV_250_RefSB_Uncert_Indexes';
    
    earthView250 = hdfread(fileName,'EV_250_RefSB');
    
    m250_scaledIntegers = bandAcrossAlong2AcrossAlongBand(earthView250);
    
    EV.radiance = scalesOffsets2Matrix(m250_scaledIntegers,radianceScales_250m,radianceOffset_250m);
    EV.reflectance = scalesOffsets2Matrix(m250_scaledIntegers,reflectanceScales_250m,reflectanceOffset_250m);
    
    
    
    
elseif strcmp(fileName(6:8),'HKM')
    
    % --------------------------------------------------------------
    % --------- FOR THIS RETRIEVAL WE ONLY NEED BANDS 1-7 ----------
    % --------------------------------------------------------------
    
    % retreive the scales and offsets for radiance and reflectance. BE
    % CAREFUL!! Each band has a unique scaling and offset. Make sure you have
    % applied the correct scaling. You can check by looking at the info
    % structure
    
    % --- Radiance scales for the first two bands ---
    radianceScales_250m = info.Vgroup.Vgroup(2).SDS(3).Attributes(6).Value; % output should be a vector with 2 entires for the first two bands
    radianceOffset_250m = info.Vgroup.Vgroup(2).SDS(3).Attributes(7).Value; % output should be a vector with 2 entires for the first two bands
    
    % --- Reflectance scales for the first two bands ---
    reflectanceScales_250m = info.Vgroup.Vgroup(2).SDS(3).Attributes(9).Value; % output should be a vector with 2 entires for the first two bands
    reflectanceOffset_250m = info.Vgroup.Vgroup(2).SDS(3).Attributes(10).Value; % output should be a vector with 2 entires for the first two bands
    
    % --- Radiance scales for bands 3-7 ---
    radianceScales_500m = info.Vgroup.Vgroup(2).SDS(1).Attributes(6).Value; % output should be a vector with 5 entries for the bands 3-7
    radianceOffsets_500m = info.Vgroup.Vgroup(2).SDS(1).Attributes(7).Value; % output should be a vector with 5 entries for the bands 3-7
    
    % --- Reflectance scales for bands 3-7 ---
    reflectanceScales_500m = info.Vgroup.Vgroup(2).SDS(1).Attributes(9).Value; % output should be a vector with 5 entries for the bands 3-7
    reflectanceOffsets_500m = info.Vgroup.Vgroup(2).SDS(1).Attributes(10).Value; % output should be a vector with 5 entries for the bands 3-7
    
    % then we are reading in calibrated Earth View data at 500m resolution,
    % including the 250m resolution bands aggregated to 500m resolution
    
    
    uncertainty = 'EV_500_RefSB_Uncert_Indexes';
    
    earthView250 = hdfread(fileName,'EV_250_Aggr500_RefSB');
    earthView500 = hdfread(fileName,'EV_500_RefSB');
    
    m250_scaledIntegers = bandAcrossAlong2AcrossAlongBand(earthView250);
    m500_scaledIntegers = bandAcrossAlong2AcrossAlongBand(earthView500);
    
    
    
    EV.radiance = cat(3, scalesOffsets2Matrix(m250_scaledIntegers,radianceScales_250m,radianceOffset_250m),...
        scalesOffsets2Matrix(m500_scaledIntegers,radianceScales_500m,radianceOffsets_500m));
    
    EV.reflectance = cat(3, scalesOffsets2Matrix(m250_scaledIntegers,reflectanceScales_250m,reflectanceOffset_250m), ...
        scalesOffsets2Matrix(m500_scaledIntegers,reflectanceScales_500m,reflectanceOffsets_500m));
    
    
    
elseif strcmp(fileName(6:8),'1KM')
    % Bands are converted to have 1km resolution
    
    % --------------------------------------------------------------
    % --------- FOR THIS RETRIEVAL WE ONLY NEED BANDS 1-7 ----------
    % --------------------------------------------------------------
    
    % retreive the scales and offsets for radiance and reflectance. BE
    % CAREFUL!! Each band has a unique scaling and offset. Make sure you have
    % applied the correct scaling. You can check by looking at the info
    % structure
    
    % --- Radiance scales for the first two bands ---
    radianceScales_250m = info.Vgroup.Vgroup(2).SDS(5).Attributes(6).Value; % output should be a vector with 2 entires for the first two bands
    radianceOffset_250m = info.Vgroup.Vgroup(2).SDS(5).Attributes(7).Value; % output should be a vector with 2 entires for the first two bands
    
    % --- Reflectance scales for the first two bands ---
    reflectanceScales_250m = info.Vgroup.Vgroup(2).SDS(5).Attributes(9).Value; % output should be a vector with 2 entires for the first two bands
    reflectanceOffset_250m = info.Vgroup.Vgroup(2).SDS(5).Attributes(10).Value; % output should be a vector with 2 entires for the first two bands
    
    % --- Radiance scales for bands 3-7 ---
    radianceScales_500m = info.Vgroup.Vgroup(2).SDS(8).Attributes(6).Value; % output should be a vector with 5 entries for the bands 3-7
    radianceOffsets_500m = info.Vgroup.Vgroup(2).SDS(8).Attributes(7).Value; % output should be a vector with 5 entries for the bands 3-7
    
    % --- Reflectance scales for bands 3-7 ---
    reflectanceScales_500m = info.Vgroup.Vgroup(2).SDS(8).Attributes(9).Value; % output should be a vector with 5 entries for the bands 3-7
    reflectanceOffsets_500m = info.Vgroup.Vgroup(2).SDS(8).Attributes(10).Value; % output should be a vector with 5 entries for the bands 3-7
    
    % --------------------------------------------------------------
    % --- NOT USING BANDS 8-36 SO WE DONT NEED EV_1KM_RefSB data ---
    % --------------------------------------------------------------
    
    
    % then we are reading in calibrated Earth View data at 1km resolution
    % including the 250m and 500m data resolution data aggregated to 1km
    % resolution
    
    uncertainty = 'EV_1KM_RefSB_Uncert_Indexes';
    
    earthView250 = hdfread(fileName,'EV_250_Aggr1km_RefSB'); % first two modis bands aggregated to 1km resolution
    
    earthView500 = hdfread(fileName,'EV_500_Aggr1km_RefSB'); % modis bands 3-7 aggregated to 1km resolution
    
    %earthView1000 = hdfread(fileName,'EV_1KM_RefSB'); % modis bands 8 - 36 at 1km resolution
    
    m250_scaledIntegers = bandAcrossAlong2AcrossAlongBand(earthView250);
    m500_scaledIntegers = bandAcrossAlong2AcrossAlongBand(earthView500);
    %m1000_scaledIntegers = bandAcrossAlong2AcrossAlongBand(earthView1000);
    
    
    EV.radiance = cat(3, scalesOffsets2Matrix(m250_scaledIntegers,radianceScales_250m,radianceOffset_250m),...
        scalesOffsets2Matrix(m500_scaledIntegers,radianceScales_500m,radianceOffsets_500m));
    
    EV.reflectance = cat(3, scalesOffsets2Matrix(m250_scaledIntegers,reflectanceScales_250m,reflectanceOffset_250m), ...
        scalesOffsets2Matrix(m500_scaledIntegers,reflectanceScales_500m,reflectanceOffsets_500m));
    
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





