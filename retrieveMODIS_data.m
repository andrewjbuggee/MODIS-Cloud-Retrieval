%% ----- Aggregate all MODIS data -----


% by Andrew J. Buggee

%%


function modis = retrieveMODIS_data(L1B_fileName,geoFileName,L2_fileName)


% extract data from hdf files
modis.EV = readMODIS_L1B_data(L1B_fileName);
[modis.sensor,modis.solar,modis.geo] = readMODIS_geolocation(geoFileName);

%Retrive the true MODIS cloud Properties

modis.cloud = readMODIS_L2_data(L2_fileName);



end