% Collect all neccessary geolocation data for calculating cloud retrievals
% from MODIS

function [sensor,solar,geo] = readMODIS_geolocation(fileName)


    info = hdfinfo(fileName);
    
    

    % load the geolocation data from the MOD03 geolocation hdf file
    % these are the lat-long positions of MODIS pixels on Earths surface
    geo.lat = hdfread(fileName,'Latitude');
    geo.long = hdfread(fileName,'Longitude');
    
    
    % Read the scale factors for the solar geometry
    solarZenith_scale = info.Vgroup.Vgroup(2).SDS(8).Attributes(4).Value;
    solarAzimuth_scale = info.Vgroup.Vgroup(2).SDS(9).Attributes(4).Value;


    % load solar position data
    solar.azimuth = hdfread(fileName,'SolarAzimuth')*solarZenith_scale; % scale factor included
    solar.zenith = hdfread(fileName,'SolarZenith')*solarAzimuth_scale; % scale factor included
    
    
    % Read the scale factors for the sensor geometry
    sensorZenith_scale = info.Vgroup.Vgroup(2).SDS(5).Attributes(4).Value;
    sensorAzimuth_scale = info.Vgroup.Vgroup(2).SDS(6).Attributes(4).Value;
    sensorRange_scale = info.Vgroup.Vgroup(2).SDS(7).Attributes(4).Value;
    
    % load satellite position data
    sensor.height = hdfread(fileName,'Height');
    sensor.range = hdfread(fileName,'Range')*sensorRange_scale; % scale factor included
    sensor.azimuth = hdfread(fileName,'SensorAzimuth')*sensorAzimuth_scale; % scale factor included
    sensor.zenith = hdfread(fileName,'SensorZenith')*sensorZenith_scale; % scale factor included


    % ---- MODIS Parameters -----


end

