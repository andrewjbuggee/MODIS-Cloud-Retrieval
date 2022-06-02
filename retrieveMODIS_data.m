%% ----- Aggregate all MODIS data -----


% by Andrew J. Buggee

%%


function [modis,L1B_fileNames] = retrieveMODIS_data(folderName)

% -----------------------------------------------------
% ----- Check to see if the folder name is valid ------
% -----------------------------------------------------


files = dir([folderName,'*.hdf']); % find all files that end in .hdf

% check to see if we found any files!
if length(files)==0
    error([newline,'There are not files in the folder provided!', newline])
end


for ii = 1:length(files)
    
    file_ii = files(ii).name;
    
    if strcmp(file_ii(1:5),'MOD02') == true || strcmp(file_ii(1:5), 'MYD02') == true
        
        if strcmp(file_ii(1:8),'MOD02HKM') == true || strcmp(file_ii(1:8),'MYD02HKM') == true
            % 500m resolution calibrated data
            
            L1B_fileNames{ii} = file_ii;
            modis.EV = readMODIS_L1B_data(file_ii);
            
        elseif strcmp(file_ii(1:8),'MOD021KM') == true || strcmp(file_ii(1:8),'MYD021KM') == true
            % 1km resolution data, which is the same resolution as the
            % cloud microphysical retrievals
            
            L1B_fileNames{ii} = file_ii;
            modis.EV = readMODIS_L1B_data(file_ii);
            
        end
        
        
        
    elseif strcmp(file_ii(1:5),'MOD03') == true || strcmp(file_ii(1:5), 'MYD03') == true
        
        % extract geolocation data from hdf files
        [modis.sensor,modis.solar,modis.geo] = readMODIS_geolocation(file_ii);
        
    elseif strcmp(file_ii(1:5),'MOD06') == true || strcmp(file_ii(1:5), 'MYD06') == true
        
        % Retrive the true MODIS cloud Properties
        
        modis.cloud = readMODIS_L2_data(file_ii);
        
    end
    
    
    
end


end