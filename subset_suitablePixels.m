%% load a subset of suitable pixels and save the indices used



% By Andrew J. Buggee
%%

function pixels2use = subset_suitablePixels(inputs,modis)

suitablePixels_fileName = [inputs.savedCalculations_folderName,'suitablePixels.mat'];

folderName2Save = inputs.savedCalculations_folderName; % where to save the indices
numPixels2Calculate = inputs.pixels.num_2calculate; % number of pixels to use in our calcualtions

load(suitablePixels_fileName,'pixels');


% for all calculations, we need the 1km resolution pixel locations. So that
% is what we will use. We only use the 500 meter resolution pixel locaiton
% to show which pixels in the EV data set we used.

numSuitablePixels = length(pixels.res1km.row);


if numSuitablePixels > numPixels2Calculate
    
    rand_indices = randi(numSuitablePixels,numPixels2Calculate,1); % generate random numbers to choose pixels
    
    pixels2use.res1km.row = pixels.res1km.row(rand_indices); % row positions
    pixels2use.res1km.col = pixels.res1km.col(rand_indices); % column positions
    
    

    pixels2use.res1km.size = pixels.res1km.size;
    %pixels2use.res500m.size = pixels.res500m.size;
    
    % lets map 1 km pixels to 500 meter pixels location

    %[pixels2use.res500m.row, pixels2use.res500m.col] = cartesian_mapping_1000_to_500_pixels(pixels2use);


    
elseif numSuitablePixels < numPixels2Calculate
    
    % if we ask for more than there are, just use them all
    
    disp(['Number of Suitable Pixels is less than what you asked for.',...
        'Writing INP files for all suitabl epixels...']);
    
    pixels2use = pixels;
    

    
    
elseif numSuitablePixels == numPixels2Calculate
    
    pixels2use = pixels;
    
    
    
    
else
    
    error('Something is wrong with the number of pixels desired')
    
end


% --- Load the geometry settings for each pixel ---

pixel_rows = pixels2use.res1km.row;
pixel_cols = pixels2use.res1km.col;

for ii = 1:length(pixel_rows)

    pixels2use.res1km.geometry.sza(ii) = modis.solar.zenith(pixel_rows(ii),pixel_cols(ii));
    pixels2use.res1km.geometry.saz(ii) = modis.solar.azimuth(pixel_rows(ii),pixel_cols(ii));
    
    % we need the cosine of the zenith viewing angle
    pixels2use.res1km.geometry.umu(ii) = round(cosd(double(modis.sensor.zenith(pixel_rows(ii),pixel_cols(ii)))),3); % values are in degrees
    pixels2use.res1km.geometry.phi(ii) = modis.sensor.azimuth(pixel_rows(ii),pixel_cols(ii));
    
end


% Save pixels2use and inputs
save([folderName2Save,inputs.saveCalculations_fileName],'pixels2use','inputs')




end