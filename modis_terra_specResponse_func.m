%% ----- read the spectral response function -----

% Response functions are stored in the rich text file (.rtf) labeled
% 'modis_spectral_response_func.rtf'

% This file contains the spectral response functions below 3 microns

% By Andrew John Buggee

%%

function spec_response = modis_terra_specResponse_func(band_number)

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
    
    

% define the filename
filename = 'modis_terra_spectral_response_func.txt';

% determine the file formatting properties
file_prop = detectImportOptions(filename);

% define the columns according to the MODIS band number
var_names = string({'Wavelength (nm)', 'Band 8', 'Band 9', 'Band 3',...
                   'Band 10', 'Band 11', 'Band 12', 'Band 4', 'Band 1',...
                   'Band 13', 'Band 14', 'Band 15', 'Band 2', 'Band 16', ...
                   'Band 5', 'Band 6', 'Band 7'}); %#ok<STRCLQT> 

% read in table
data = readtable(filename, file_prop);

% reset the variable names
data.Properties.VariableNames = var_names;


% The last column needs to be fixed. It's read in as a cell vector of
% string characters



% % open the file for reading
% file_id = fopen(filename, 'r');   % 'r' tells the function to open the file for reading
% 
% format_spec = '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f';                                  % two floating point numbers
% data = textscan(file_id, format_spec, 'Delimiter',' ',...
% 'MultipleDelimsAsOne',1, 'CommentStyle','/', 'EndOfLine','\');

%% Read in the correct response function

% Store the wavelength data in a vector. Values associated with non-zero
% responses will be thrown out

wavelength = data.("Wavelength (nm)");

% define an empty cell aray 
spec_response = cell(1, length(band_number));


for nn = 1:length(band_number)

    for bb = 1:(length(data.Properties.VariableNames)-1)

        if strcmp(['Band ',num2str(band_number(nn))], data.Properties.VariableNames{bb+1})==true

            % Select this spectral response function
            data2keep = data{:,bb+1};

            % only keep the non-zero values
            index_nonZero = find(data2keep);

            spec_response{nn}(:,1) = wavelength(index_nonZero);
            spec_response{nn}(:,2) = data2keep(index_nonZero);


        end
    end

end


end



