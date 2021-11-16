%% ----- Extract MODIS calculated values as truth -----



% By Andrew J. Buggee

%%

function [truthEstimateTable, reflectanceTable] = gatherTruthEstimateVals(modis, minVals, inputs)

% extract inputs


bands2search = inputs.bands2search;
pixel_row = inputs.pixel_row;
pixel_col = inputs.pixel_col;

truthEstimateTable = table;
reflectanceTable = table;

% extract the modis computed values for reflectance to compare with my
% internal computation of reflectance
for ii = 1:length(bands2search)
    if bands2search(ii)<=2
        
        reflectanceTable.modisRefl(ii) = modis.EV.m250.reflectance(pixel_row,pixel_col,bands2search(ii));
        
    elseif bands2search(ii)>2
        
         reflectanceTable.modisRefl(ii) = modis.EV.m500.reflectance(pixel_row,pixel_col,bands2search(ii)-2);
    end
    
end


% extract the values MODIS calculates using its own algorithm
truthEstimateTable.modisR17 = modis.cloud.effRadius17(pixel_row,pixel_col);
truthEstimateTable.modisT17 = modis.cloud.optThickness17(pixel_row,pixel_col);

% gather the estimated values in the table

truthEstimateTable.estR17 = minVals.minR;
truthEstimateTable.estT17 = minVals.minT;


% compute the squared difference
truthEstimateTable.absDiffR = abs(minVals.minR - truthEstimateTable.modisR17);
truthEstimateTable.absDiffT = abs(minVals.minT - truthEstimateTable.modisT17);







end
