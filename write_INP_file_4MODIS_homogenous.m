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

%% DEFINE UVSPEC INP FILE INPUTS


% ---------------------------------------------------
% -------- Define the solar flux file to use --------
% ---------------------------------------------------
if inputs.solarFlux_resolution==1
    % use 1nm resolution
    solarFlux_file = 'kurudz_1.0nm.dat';

elseif inputs.solarFlux_resolution==0.1
    % use 0.1nm resoltion file
    solarFlux_file = 'kurudz_0.1nm.dat';

else
    error([newline, 'I dont recognize the solar flux file. Only options are the 1nm and 0.1nm Kuruz files',newline])
    
end



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


% define the MODIS bands to run using the spectral response functions
lambda = zeros(length(inputs.bands2run), 2);

% now, we need to may need to interpolate these curves, depending on the
% resolution of the solar flux file that we chose
% place it on the same wavelength grid as the solar file setting
% do this using interpolation
if inputs.solarFlux_resolution==0.1

    % interpolate over all bands that are being run
    for bb = 1:length(inputs.bands2run)
        % define the new wavelength vector
        new_wl_vec = modis.spec_response{bb}(1,1):0.1:modis.spec_response{bb}(end,1);

        interp_spec_response = interp1(modis.spec_response{bb}(:,1), modis.spec_response{bb}(:,2),new_wl_vec);

        % rewrite the spec_response cell array
        modis.spec_response{bb} = [new_wl_vec', interp_spec_response'];

        % Now we can define the wavelength bounds for each band
        % The wavelength grid and step size is defined by the solar file
        lambda(bb,:) = [new_wl_vec(1), new_wl_vec(end)];      % nm - midpoint and boundaries for each MODIS band we wish to model

    end

elseif inputs.solarFlux_resolution==1

    % in this case, the solar file jumps every nm and we don't need to
    % interpoalte the spectral response functions
    for bb = 1:length(inputs.bands2run)

        lambda(bb,:) = [modis.spec_response{bb}(1,1), modis.spec_response{bb}(end,1)];      % nm - midpoint and boundaries for each MODIS band we wish to model

    end

end




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

% For now, lets only use the midpoint of MODIS band 1
lambda_forTau = round((lambda(1,2) - lambda(1,1))/2 + lambda(1,1));                % nm

% define the day of the year and correct the Earth sun distance
day_of_year = inputs.day_of_year;


% ------------------------------------------------------
% ---------- Lets define the % of cloud cover ----------
% ------------------------------------------------------
%cloud_cover = 0.775;       % this value seemed to work
cloud_cover = 1;
% ------------------------------------------------------


% we can create all the water files ahead of the INP file creation. Then we
% will step through each component and create an INP file

% ------------------------------------------------------
% --- WE COULD USE CLOUD HEIGHT DERIVED FROM MODIS?! ---
% ------------------------------------------------------

% cloudTopHeight = modis.cloud.topHeight(pixels2use.res1km.linearIndex)./1e3;         % km - MODIS estiamte cloud top height  
% z_topBottom = [cloudTopHeight, cloudTopHeight-0.5];                                 % hard code all clouds to be 500 meters thick
% % Find the locations where clouds are at the surface
% index0 = z_topBottom(:,2)<0;
% z_topBottom(index0,2) =0;
% z_topBottom = z_topBottom';                                             % needs to be in an array with 2 rows
% ---------------------------------------
% For now lets hard code the cloud height
% ---------------------------------------
% cloud top height and thickness is set to be the same values 
% described in: "Overview of the MODIS Collection 6 Cloud Optical 
% Property (MOD06) Retrieval Look-up Tables" Amarasinghe
% et. al 2017
z_topBottom = [9, 8];                     % km - altitude above surface for the cloud top and cloud bottom
% ---------------------------------------
% ---------------------------------------

% MODIS integrates over a modified gamma droplet distribution
dist_str = inputs.clouds.distribution_type;
% libRadTran stats in their manual that a distribution variance value of 7
% is adequate for liquid water clouds. The distribution variance needs to
% be the same size as the effective radius vector
dist_var = repmat(inputs.clouds.distribution_variance,size(re_mat,1), size(re_mat,2));

% Since this is writing files to do the TBLUT retrieval, 
% we assume clouds are vertical homogeneous
vert_homogeneity_str = 'vert-homogeneous';
parameterization_str = inputs.clouds.wc_parameterization;        % Tells the function how to compute LWC

wc_filename = write_wc_file(re_mat, tau_mat(:), z_topBottom, lambda_forTau, dist_str, dist_var,...
    vert_homogeneity_str, parameterization_str);

% rehape so we can step through values for r and tau
wc_filename = reshape(wc_filename, length(re),length(tau_c));

% define the parameterization used to convert water cloud properties to
% optical properties
wc_parameterization = inputs.clouds.wc_properties;

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
    
    % ----- Define the Solar Geometry -------
    sza = modis.solar.zenith(pixel_row(pp),pixel_col(pp));
    phi0 = 180 + modis.solar.azimuth(pixel_row(pp),pixel_col(pp));      % degree - this is how we map MODIS azimuth of the sun to the LibRadTran measurement
    

    % ------ Define the Viewing Geometry ------
    % we need the cosine of the zenith viewing angle
    % positive values solve for upwelling radiance, where the sensor is
    % defined to be looking down towrads the Earth's surface. negative
    % values solve for downwelling radiance, where the sensor is looking
    % upwards towards the sky
    umu = round(cosd(double(modis.sensor.zenith(pixel_row(pp),pixel_col(pp)))),4); % values are in degrees

    % to properly map the azimuth angle onto the reference plane used by
    % libRadTran, we need an if statement
    if modis.sensor.azimuth(r,c)<0
        phi = 360 + modis.sensor.azimuth(pixel_row(pp),pixel_col(pp));
    else
        phi = modis.sensor.azimuth(pixel_row(pp),pixel_col(pp));
    end
    
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
    
    
    % For each pixel, use the MODIS retrieved column water vapor and
    % overide the US standard atmosphere profile
    column_vapor = modis.vapor.col_nir(pixel_row(pp),pixel_col(pp))*10;         % mm - column of total precipitable water over a square mm
    
    
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
                
                
                % Define the location and filename of the US standard
                % atmosphere that will be used for this analysis
                formatSpec = '%s %s %5s %s \n';
                fprintf(fileID, formatSpec,'atmosphere_file','../data/atmmod/afglus.dat',' ', '# Location of atmospheric profile to use');
                
                % Define the location and filename of the extraterrestrial solar source
                formatSpec = '%s %s %5s %s \n\n';
                fprintf(fileID, formatSpec,'source solar',['../data/solar_flux/',solarFlux_file], ' ', '# Bounds between 250 and 10000 nm');
                
                % Define the location and filename of the extraterrestrial solar source
                formatSpec = '%s %f %5s %s \n\n';
                fprintf(fileID, formatSpec,'day_of_year', day_of_year, ' ', '# accounts for changing Earth-Sun distance');

                % Define the ozone column
%                 formatSpec = '%s %s %s %s %5s %s \n';
%                 fprintf(fileID, formatSpec,'mol_modify','O3', '300.','DU', ' ', '# Set ozone column');

                % Define the water vapor column
                formatSpec = '%s %s %s %s %5s %s \n';
                fprintf(fileID, formatSpec,'mol_modify','H2O', num2str(column_vapor), 'MM', ' ', '# Total Precipitable Water');
                
                % Define the surface albedo
                formatSpec = '%s %s %5s %s \n\n';
                fprintf(fileID, formatSpec,'albedo','0.0600', ' ', '# Surface albedo of the ocean');
                
                
                % Define the water cloud file
                formatSpec = '%s %s %5s %s \n';
                fprintf(fileID, formatSpec,'wc_file 1D', ['../data/wc/',wc_filename{rr,tt}], ' ', '# Location of water cloud file');
                
                
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
                fprintf(fileID, formatSpec,'wavelength', lambda(bb,1), lambda(bb,2), ' ', '# Wavelength range');
                
                
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

