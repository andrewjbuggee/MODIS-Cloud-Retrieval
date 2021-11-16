%% load a subset of suitable pixels and save the indices used



% By Andrew J. Buggee
%%

function pixels2use = load_suitablePixels(inputs)

suitablePixels_fileName = [inputs.savedCalculations_folderName,'suitablePixels.mat'];

folderName2Save = inputs.savedCalculations_folderName; % where to save the indices
numPixels = inputs.pixels.num_INP_files_2write; % number of pixels to use in our calcualtions

load(suitablePixels_fileName,'pixels');


% for all calculations, we need the 1km resolution pixel locations. So that
% is what we will use. We only use the 500 meter resolution pixel locaiton
% to show which pixels in the EV data set we used.

numSuitablePixels = length(pixels.res1km.row);


if numSuitablePixels > numPixels
    
    rand_indices = randi(numSuitablePixels,numPixels,1); % generate random numbers to choose pixels
    
    pixels2use.res1km.row = pixels.res1km.row(rand_indices); % row positions
    pixels2use.res1km.col = pixels.res1km.col(rand_indices); % column positions
    
    

    pixels2use.res1km.size = pixels.res1km.size;
    pixels2use.res500m.size = pixels.res500m.size;
    
    % lets map 1 km pixels to 500 meter pixels location

    [pixels.res500m.row, pixels.res500m.col] = cartesian_mapping_1000_to_500_pixels(pixels);

    save([folderName2Save,inputs.saveCalculations_fileName],'pixels2use','inputs')
    
elseif numSuitablePixels < numPixels
    
    % if we ask for more than there are, just use them all
    
    disp(['Number of Suitable Pixels is less than what you asked for.',...
        'Writing INP files for all suitabl epixels...']);
    
    pixels2use = pixels;
    
    save([folderName2Save,inputs.saveCalculations_fileName],'pixels2use','inputs')
    
    
elseif numSuitablePixels == numPixels
    
    pixels2use = pixels;
    
    save([folderName2Save,inputs.saveCalculations_fileName],'pixels2use','inputs')
    
    
else
    
    error('Something is wrong with the number of pixels desired')
    
end






end