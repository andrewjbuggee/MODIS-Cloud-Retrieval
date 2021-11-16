% Collect all neccessary geolocation data for calculating cloud retrievals
% from MODIS

function [sensor,solar,geo] = readMODIS_geolocation(fileName)




    % load the geolocation data from the MOD03 geolocation hdf file
    % these are the lat-long positions of MODIS pixels on Earths surface
    geo.lat = hdfread(fileName,'Latitude');
    geo.long = hdfread(fileName,'Longitude');
    

    % load solar position data
    solar.azimuth = hdfread(fileName,'SolarAzimuth')./100; % scale factor included
    solar.zenith = hdfread(fileName,'SolarZenith')./100; % scale factor included
    
    % load satellite position data
    sensor.height = hdfread(fileName,'Height');
    sensor.range = hdfread(fileName,'Range').*25; % scale factor included
    sensor.azimuth = hdfread(fileName,'SensorAzimuth')./100; % scale factor included
    sensor.zenith = hdfread(fileName,'SensorZenith')./100; % scale factor included


    % ---- MODIS Parameters -----


end

