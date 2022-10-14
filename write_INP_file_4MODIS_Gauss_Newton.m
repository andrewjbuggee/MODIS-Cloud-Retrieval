%% This function will write a .INP uvspec file for LibRadTran



%%

function inpNames = write_INP_file_4MODIS_Gauss_Newton(GN_inputs, modisInputs, pixel_row, pixel_col, modis, wc_filename)

% ------------------------------------------------
% ---------- INPUTS AND FUNCTION SET UP ----------
% ------------------------------------------------


% what computer are we using?

% a template file has been set up to be edited and saved as a new file
% determine which computer is being used
userName = whatComputer;

if strcmp(userName,'anbu8374')
    
    libRadTran_path = ['/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4'];
    
elseif strcmp(userName,'andrewbuggee')
    
    libRadTran_path = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/',...
        'LibRadTran/libRadtran-2.0.4'];
    
else
    error('I dont recognize this computer user name')
end


addpath(libRadTran_path);

%%

% for each spectral bin, we have an image on the ground composed of 2030 *
% 1354 pixels. The swath on the ground is large enough that the solar
% zenith and solar azimuth change significantly. Ideally at each spectral bin, and
% each pixel, we could calculate the reflectance function and how it varies
% with changing tau and change re.


newFolder = [libRadTran_path,'/',modisInputs.INP_folderName]; % where the newly created .inp files will be saved

% If the folder doesn't exist, create it
if ~exist(newFolder, 'dir')
    mkdir(newFolder)
end

% -------------------------------------------------
% -------------------------------------------------
% define the MODIS bands of interest

lambda = modisBands(GN_inputs.bands2use);      % nm - midpoint and boundaries for each MODIS band we wish to model
% -------------------------------------------------
% -------------------------------------------------


% define the parameterization used to convert water cloud properties to
% optical properties
wc_parameterization = modisInputs.flags.wc_properties;

% we have 4 variables that can change for each INP file

%   1) pixel
%   2) modis band
%   3) re
%   4) tau_c

% ------------------------------------------------------
% ---------- Lets define the % of cloud cover ----------
% ------------------------------------------------------

if strcmp(modisInputs.modisDataFolder(96:end), '/2008_11_11_1850/')==true

    % For 11/11/2008 - 14:30 data
    cloud_cover = 0.795;

elseif strcmp(modisInputs.modisDataFolder(96:end), '/2008_11_11_1430/')==true

    % For 11/11/2008 - 14:30 data
    cloud_cover = 0.725;

elseif strcmp(modisInputs.modisDataFolder(96:end), '/2008_11_09/')==true

    % For 11/09/2008 data
    cloud_cover = 0.70;

else

    cloud_cover = 0.775;

end








% ------------------------------------------------------



% write INP files for all 7 MODIS bands. This function writes files for a
% single pixel at a time, for now...
numPixels = length(pixel_row);
numBands = length(GN_inputs.bands2use);

inpNames = cell(numPixels, numBands);

for pp = 1:numPixels
    
    sza = modis.solar.zenith(pixel_row(pp),pixel_col(pp));
    phi0 = modis.solar.azimuth(pixel_row(pp),pixel_col(pp));
    
    % we need the cosine of the zenith viewing angle
    umu = round(cosd(double(modis.sensor.zenith(pixel_row(pp),pixel_col(pp)))),3); % values are in degrees
    phi = modis.sensor.azimuth(pixel_row(pp),pixel_col(pp));
    
    % Modis defines the azimuth viewing angle as [0,180]
    % and [-180,0], whereas libradtran defines the azimuth
    % angle as [0,360]. So we need to make this adjustment
    
    if phi<0
        phi = phi+360;
    end
    
    % Modis defines the azimuth viewing angle as [0,180]
    % and [-180,0], whereas libradtran defines the azimuth
    % angle as [0,360]. So we need to make this adjustment
    
    if phi0<0
        phi0 = phi0+360;
    end
    
    % create the begining of the file name string
    fileBegin = ['GN_pixel_',num2str(pixel_row(pp)),'r_',num2str(pixel_col(pp)),'c_sza_',num2str(sza),'_saz_',num2str(phi0),'_band_'];
    
    for bb = 1:numBands
        
        band_num = GN_inputs.bands2use(bb);        % modis band number that defines the upper and lower wavelength boundaries used to compute the equation of radiative transfer
        
        
        
        
        % redefine the old file each time
        inpNames{pp,bb} = [fileBegin,num2str(band_num),'.INP'];
        
        % ------------------------------------------------------------
        % --------------------- WRITE INP FILE -----------------------
        % ------------------------------------------------------------
        
        
        
        % Create the water cloud file
        fileID = fopen([newFolder,inpNames{pp,bb}], 'w');
        
        % fprintf writes lines in our text file from top to botom
        % wc.DAT files are written with the higher altitudes at the top, and the
        % surface at the bottom
        
        % to write column vectors in a text file, we have to store them as row
        % vectors
        
        % The first argument is the format specification
        % The designates what type of characters to print, such as floating point
        % arithmetic or string. The number before the character designates the
        % MINIMUM number of characters to print
        
        
        % Define which RTE solver to use
        formatSpec = '%s %s %5s %s \n';
        fprintf(fileID, formatSpec,'rte_solver','disort',' ', '# Radiative transfer equation solver');
        
        
        % Define the number of streams to keep track of when solving the equation
        % of radiative transfer
        formatSpec = '%s %s %5s %s \n\n';
        fprintf(fileID, formatSpec,'number_of_streams','6',' ', '# Number of streams');
        
        
        % Define the location and filename of the atmopsheric profile to use
        formatSpec = '%s %s %5s %s \n';
        fprintf(fileID, formatSpec,'atmosphere_file','../data/atmmod/afglus.dat',' ', '# Location of atmospheric profile to use');
        
        % Define the location and filename of the extraterrestrial solar source
        formatSpec = '%s %s %5s %s \n\n';
        fprintf(fileID, formatSpec,'source solar','../data/solar_flux/kurudz_1.0nm.dat', ' ', '# Bounds between 250 and 10000 nm');
        
        % Define the ozone column
        formatSpec = '%s %s %s %s %5s %s \n';
        fprintf(fileID, formatSpec,'mol_modify','O3', '300.','DU', ' ', '# Set ozone column');
        
        % Define the surface albedo
        formatSpec = '%s %s %5s %s \n\n';
        fprintf(fileID, formatSpec,'albedo','0.0600', ' ', '# Surface albedo of the ocean');
        
        
        % Define the water cloud file
        formatSpec = '%s %s %5s %s \n';
        fprintf(fileID, formatSpec,'wc_file 1D', ['../data/wc/',wc_filename{1}], ' ', '# Location of water cloud file');
        
        
        % Define the percentage of horizontal cloud cover
        % This is a number between 0 and 1
        formatSpec = '%s %s %5s %s \n';
        fprintf(fileID, formatSpec,'cloudcover wc', num2str(cloud_cover), ' ', '# Cloud cover percentage');
        
        % Define the technique or parameterization used to convert liquid cloud
        % properties of r_eff and LWC to optical depth
        formatSpec = '%s %s %5s %s \n\n';
        fprintf(fileID, formatSpec,'wc_properties', wc_parameterization, ' ', '# optical properties parameterization technique');
        
        % Define the wavelengths for which the equation of radiative transfer will
        % be solve
        formatSpec = '%s %f %f %5s %s \n\n';
        fprintf(fileID, formatSpec,'wavelength', lambda(bb,2), lambda(bb,3), ' ', '# Wavelength range');
        
        
        % Define the sensor altitude
        formatSpec = '%s %s %5s %s \n';
        fprintf(fileID, formatSpec,'zout', 'toa', ' ', '# Sensor Altitude');
        
        % Define the solar zenith angle
        formatSpec = '%s %f %5s %s \n';
        fprintf(fileID, formatSpec,'sza', sza, ' ', '# Solar zenith angle');
        
        % Define the solar azimuth angle
        formatSpec = '%s %f %5s %s \n';
        fprintf(fileID, formatSpec,'phi0', phi0, ' ', '# Solar azimuth angle');
        
        % Define the cosine of the zenith viewing angle
        formatSpec = '%s %f %5s %s \n';
        fprintf(fileID, formatSpec,'umu', umu, ' ', '# Cosine of the zenith viewing angle');
        
        % Define the azimuth viewing angle
        formatSpec = '%s %f %5s %s \n\n';
        fprintf(fileID, formatSpec,'phi', phi, ' ', '# Azimuthal viewing angle');
        
        
        % Set the error message to quiet of verbose
        formatSpec = '%s';
        fprintf(fileID, formatSpec,'quiet');
        
        
        fclose(fileID);   
        
        
    end
    
    
    
    
end



end

