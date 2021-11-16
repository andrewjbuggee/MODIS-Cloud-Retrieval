% ----- Determine the MODIS Bands -----

function bandVals = modisBandsCenter(band_number)

% Check inputs 

    if isnumeric(band_number)
        
    else
        error('input must be a numeric entry')
        
    end
    
    if length(band_number)>36
        error('MODIS has 36 spectral bands. You requested more than 36. You may pull any of of them')
    end
    
    if sum(band_number>36)>0
        error('You requested a band number that is higher than 36, which doesnt exist')
    end
    
    

filename = 'modis_bands.txt';

bandData = readtable(filename);


lowVals = bandData.Var2(band_number);
upperVals = bandData.Var3(band_number);

midVals = (upperVals - lowVals)/2 + lowVals;

bandVals = midVals;

end



