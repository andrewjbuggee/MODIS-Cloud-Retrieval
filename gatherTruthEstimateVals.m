%% ----- Extract MODIS calculated values as truth -----



% By Andrew J. Buggee

%%

function [truthTable, reflectanceTable] = gatherTruthEstimateVals(modis, minVals, inputs,pixels2use)

% extract inputs

bands2run = inputs.bands2run;
bands2plot = inputs.bands2plot;
num_pixels = inputs.pixels.num_2calculate;
bands2search = inputs.bands2search;

% 1km resolution pixels
pixel_row_1 = pixels2use.res1km.row;
pixel_col_1 = pixels2use.res1km.col;

% 500 meter resolution pixels
pixel_row_500 = pixels2use.res500m.row;
pixel_col_500 = pixels2use.res500m.col;

% save calculations
saveCalcs_filename = inputs.saveCalculations_fileName;


% create reflectance function table
reflectanceTable = table;
truthTable = table;


% create cell array for band names
Band = cell(length(bands2search),1);


% extract the modis computed values for reflectance to compare with my
% internal computation of reflectance

% SWITCH TO 1KM DATA
% but for now, lets propose a simple algebraic fix so taht we select
% the proper pixel in our 500 meter data set
% lets take an average of the four pixels that map to the 500 meter
% data
pixelIndex_500m = 1:4:((num_pixels-1)*4 +1);

for pp = 1:num_pixels
    
    
    for ii = 1:length(bands2search)
        
        Band{ii} = [num2str(modisBandsCenter(bands2search(ii))),' nm'];
        
        if bands2search(ii)<=2
            
            % lets grab the 4 pixels that the 1 km data maps to in the 500
            % meter data space. Then take the average value
            four_vals_index = pixelIndex_500m(pp):1:(pixelIndex_500m(pp)+3);
            for kk = 1:length(four_vals_index)
                
                hold4_values(kk) = modis.EV.m250.reflectance(pixel_row_500(four_vals_index(kk)),pixel_col_500(four_vals_index(kk)),bands2search(ii));
            end
            
            reflectanceTable.modisRefl(ii,pp) = mean(hold4_values);
            
        elseif bands2search(ii)>2
            
            % lets grab the 4 pixels that the 1 km data maps to in the 500
            % meter data space. Then take the average value
            four_vals_index = pixelIndex_500m(pp):1:(pixelIndex_500m(pp)+3);
            for kk = 1:length(four_vals_index)
                
                hold4_values(kk) = modis.EV.m500.reflectance(pixel_row_500(four_vals_index(kk)),pixel_col_500(four_vals_index(kk)),bands2search(ii)-2);
            end
            reflectanceTable.modisRefl(ii,pp) = mean(hold4_values);
            
            
        end
        
    end
    
    
    
    % extract the values MODIS calculates using its own algorithm
    truthTable.modisR17(pp) = modis.cloud.effRadius17(pixel_row_1(pp),pixel_col_1(pp));
    truthTable.modisT17(pp) = modis.cloud.optThickness17(pixel_row_1(pp),pixel_col_1(pp));
    
    % gather the estimated values in the table
    
    truthTable.estR17(pp) = minVals.minR(pp);
    truthTable.estT17(pp) = minVals.minT(pp);
    
    % compute the absolute difference
    truthTable.absDiffR(pp) = abs(minVals.minR(pp) - truthTable.modisR17(pp));
    truthTable.absDiffT(pp) = abs(minVals.minT(pp) - truthTable.modisT17(pp));
    
    
    % compute the percent difference
    truthTable.percentDiffR(pp) = abs(1 - minVals.minR(pp)/truthTable.modisR17(pp)) * 100;
    truthTable.percentDiffT(pp) = abs(1 - minVals.minT(pp)/truthTable.modisT17(pp)) * 100;
end

reflectanceTable.Bands = Band;




save(saveCalcs_filename,"truthTable","reflectanceTable",'-append'); % save inputSettings to the same folder as the input and output file






end
