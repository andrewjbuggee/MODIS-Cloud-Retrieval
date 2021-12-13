%% ----- Aggregate all MODIS data -----


% by Andrew J. Buggee

%%


function [modis,L1B_500m_fileName] = retrieveMODIS_data(folderName)


files = dir([folderName,'*.hdf']); % find all files that end in .hdf

for ii = 1:length(files)
    
    file_ii = files(ii).name;
    
    if strcmp(file_ii(1:5),'MOD02') == true
        
        if strcmp(file_ii(1:8),'MOD02HKM') == true
            
            L1B_500m_fileName = file_ii;
            modis.EV = readMODIS_L1B_data(file_ii);
            
        elseif strcmp(file_ii(1:8),'MOD021KM') == true
            
            % RIGHT NOW, DO NOTHING
        end
        
    elseif strcmp(file_ii(1:5),'MOD03') == true
        
        % extract geolocation data from hdf files
        [modis.sensor,modis.solar,modis.geo] = readMODIS_geolocation(file_ii);
        
    elseif strcmp(file_ii(1:5),'MOD06') == true
        
        %Retrive the true MODIS cloud Properties
        
        modis.cloud = readMODIS_L2_data(file_ii);
        
    end
    
    
    
end


end