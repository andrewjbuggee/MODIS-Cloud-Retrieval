% grab modis reflectacne data for the bands of interest


% By Andrew J. Buggee
%%

function modisRefl = grab_modis_reflectance(modis,inputs)

bands2run = inputs.bands2run;

modisRefl = [];

for ii = 1:length(bands2run)
    
    modisRefl = cat(3,modisRefl,modis.EV1km.reflectance(:,:,bands2run(ii)));
    
    
end




end