%% This function will write a .INP uvspec file for LibRadTran



%%

function inpNames = write_INP_file_4MODIS_homogenous(inputs, pixels2use, modis)

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


pixel_row = pixels2use.res1km.row; % for everything I need in this code, we use the 1km resolution pixel locations
pixel_col = pixels2use.res1km.col; %
newFolder = [libRadTran_path,'/',inputs.INP_folderName]; % where the newly created .inp files will be saved

% If the folder doesn't exist, create it
if ~exist(newFolder, 'dir')
       mkdir(newFolder)
end


% define the MODIS bands of interest
lambda = modisBands(inputs.bands2run);      % nm - midpoint and boundaries for each MODIS band we wish to model


% MODIS only considers homogenous plane parallel clouds. Lets construct the
% re matrix needed to create homogenous water clouds using write_wc_file
re = inputs.re;                 % microns - values of re that we wish to model
tau_c = inputs.tau_c;           % values of tau that we wish to model

re_mat = repmat(inputs.re,1,length(inputs.tau_c));      % homogenous re matrix
tau_mat = repmat(inputs.tau_c,length(inputs.re),1);  % we need to run every re at every value of tau

% % the above matrices must be rerun for every modis band!
% re = repmat(re,1,length(inputs.bands2run));
% tau_c = repmat(tau_c, 1, length(inputs.bands2run));

% create the lambda matrix with all values needed to create INP files

% for the write_wc_file function we need to define the wavelength that is
% used to compute the optical depth. It should be bands 1 3 or 4.

% For now, lets only use MODIS band 1
lambda_forTau = lambda(1,1);                % nm

% we can create all the water files ahead of the INP file creation. Then we
% will step through each component and create an INP file

% ------------------------------------------------------
% --- WE COULD USE CLOUD HEIGHT DERIVED FROM MODIS?! ---
% ------------------------------------------------------

% For now lets hard code the cloud height
z_topBottom = [1.5, 1];         % km - altitude above surface for the cloud top and cloud bottom

% lets only model the monodispersed clouds
dist_str = 'mono';
homogeneity_str = 'homogeneous';

wc_filename = write_wc_file(re_mat, tau_mat(:), z_topBottom, lambda_forTau, dist_str, homogeneity_str);

% rehape so we can step through values for r and tau
wc_filename = reshape(wc_filename, length(re),length(tau_c));

% define the parameterization used to convert water cloud properties to
% optical properties
wc_parameterization = inputs.flags.wc_properties;

% we have 4 variables that can change for each INP file

%   1) pixel
%   2) modis band
%   3) re
%   4) tau_c



% step through each pixel, and feed all re values, tau values, and modis
% bands of interest


% step through each band, each effective raidus and each optical depth
inpNames = cell(length(pixel_row), length(re),length(tau_c),length(inputs.bands2run));

for pp = 1:length(pixel_row)
    
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
    fileBegin = ['pixel_',num2str(pixel_row(pp)),'r_',num2str(pixel_col(pp)),'c_sza_',num2str(sza),'_saz_',num2str(phi0),'_band_'];
    
    for bb = 1:length(inputs.bands2run)
        
        band_num = inputs.bands2run(bb);        % modis band number that defines the upper and lower wavelength boundaries used to compute the equation of radiative transfer
        
        
        for rr = 1:length(inputs.re)
            
            
            
            for tt = 1:length(inputs.tau_c)
                
                
                % redefine the old file each time
                inpNames{pp,rr,tt,bb} = [fileBegin,num2str(band_num),'_r_',num2str(re(rr)),'_T_',num2str(tau_c(tt)),'.INP'];
                
                % ------------------------------------------------------------
                % --------------------- WRITE INP FILE -----------------------
                % ------------------------------------------------------------
                
                
                
                % Create the water cloud file
                fileID = fopen([newFolder,inpNames{pp,rr,tt,bb}], 'w');
                
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
                fprintf(fileID, formatSpec,'wc_file 1D', ['../data/wc/',wc_filename{rr,tt}], ' ', '# Location of water cloud file');
                
                
                % Define the percentage of horizontal cloud cover
                % This is a number between 0 and 1
                formatSpec = '%s %s %5s %s \n';
                fprintf(fileID, formatSpec,'cloudcover wc', '1.0', ' ', '# Cloud cover percentage');
                
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
    
    
    
    
    
    
    
end



end

