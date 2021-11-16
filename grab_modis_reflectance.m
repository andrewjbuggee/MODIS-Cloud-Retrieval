% grab modis reflectacne data for the bands of interest


% By Andrew J. Buggee
%%

function modisR = grab_modis_reflectance(modis,inputs)

bands2run = inputs.bands2run;

modisR = [];

for ii = 1:length(bands2run)
    
    
    if bands2run(ii)<=2
        
        index = bands2run(ii);
        modisR = cat(3,modisR,modis.EV.m250.reflectance(:,:,index));
        
    elseif bands2run(ii)>2
        
        index = bands2run(ii) - 2;
        modisR = cat(3,modisR,modis.EV.m500.reflectance(:,:,index));
        
    end
    
    
end




end