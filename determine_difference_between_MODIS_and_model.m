%% Determine the avg difference between my estimate of reflectance the MODIS measurement

%By Andrew John Buggee

%% Loop through 100 different MODIS Pixels

clear variables

% Load MODIS data set

modisFolder = '/Users/anbu8374/Documents/MATLAB/HyperSpectral_Cloud_Retrieval/MODIS_Cloud_Retrieval/MODIS_data/2023_04_13/';

[modis,L1B_500m_fileName] = retrieveMODIS_data(modisFolder);

% Grab n random pixels from the suitablePixels mat file
load('/Users/anbu8374/Documents/MATLAB/HyperSpectral_Cloud_Retrieval/MODIS_Cloud_Retrieval/MODIS_data/2023_04_13/suitablePixels.mat', 'pixels')

n_pixels = 300;

% --- We only want pixles with a droplet size less than 25 ---
index_25 = modis.cloud.effRadius17(pixels.res1km.index)<25;
% --- Let's only look at pixels within a narrow reflectance range ---
% modis_reflectance = modis.EV1km.reflectance(:, :, 1);
% index_refl = modis_reflectance(pixels.res1km.index)>=0.24 & modis_reflectance(pixels.res1km.index)<=0.26;
% define the sample population with the indexes just found
population = pixels.res1km.index(index_25);
% sample from the indexes above
idx = randsample(population, n_pixels, false);

% determine the rows and columns
[row, col] = ind2sub(size(modis.cloud.effRadius17), idx);

% clear the pixels strucutre, we don't need it anymore
clear('pixels');
%% Define the parameters of the INP file


% Define the number of streams to use in your radiative transfer model
num_streams = 16;
% ------------------------------------------------------------------------


% define the wavelength range. If monochromatic, enter the same number
% twice
band1 = modisBands(1);
wavelength = [band1(2), band1(3)];              % nm - monochromatic wavelength calcualtion


% ------------------------------------------------------------------------
% --- Do you want to use the Nakajima and Tanka radiance correction? -----
use_nakajima_phaseCorrection = true;
% ------------------------------------------------------------------------

% define the atmospheric data file
atm_file = 'afglus.dat';

% define the surface albedo
albedo = 0.05;

% day of the year
day_of_year = str2double(L1B_500m_fileName{1}(15:17));

% ------------------------------------------------------------------------
% ------ Do you want to use the MODIS cloud top height estimate? ---------
use_MODIS_cloudTopHeight = true;
% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% -------------- Do you want a cloud in your model? ----------------------
yesCloud = true;

% ---- Do you want a linear adjustment to the cloud pixel fraction? ------
linear_cloudFraction = true;
% if false, define the cloud cover percentage
cloud_cover = 1;
% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% ---------- Do you want use your custom mie calculation file? -----------
use_custom_mie_calcs = false;
% ------------------------------------------------------------------------


% define the type of droplet distribution
distribution_str = 'gamma';
% define the distribution varaince
distribution_var = 10;
% define whether this is a vertically homogenous cloud or not
vert_homogeneous_str = 'vert-homogeneous';
% define how liquid water content will be computed
parameterization_str = 'mie';

% define the wavelength used for the optical depth
lambda_forTau = (wavelength(2) - wavelength(1))/2 + wavelength(1);            % nm



% Define the parameterization scheme used to comptue the optical quantities
if use_custom_mie_calcs==false
    wc_parameterization = 'mie interpolate';
else
    %wc_parameterization = '../data/wc/mie/wc.mie_test.cdf interpolate';
    wc_parameterization = '../data/wc/mie/wc.mie_test2_more_nmom.cdf interpolate';
end

% --------------------------------------------------------------
% --------------------------------------------------------------



% --------------------------------------------------------------
% --- Do you want to use the Cox-Munk Ocean Surface Model? -----
use_coxMunk = true;
wind_speed = 3;             % m/s
% --------------------------------------------------------------


% ------------------------------------------------------------------------
% --------- Do you want boundary layer aerosols in your model? -----------
yesAerosols = true;

aerosol_type = 4;               % 4 = maritime aerosols
aerosol_opticalDepth = 0.1;     % MODIS algorithm always set to 0.1
% ------------------------------------------------------------------------







%% Write each INP file using various MODIS values

% Define the folder path where all .INP files will be saved
folder2save = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/solving_modeling_discrepancy_2/';

inputName = cell(1, length(idx));
outputName = cell(1, length(idx));
wc_filename = cell(1, length(idx));

% save the modis droplet size and optical depth
re = zeros(1, length(idx));
tau_c = zeros(1, length(idx));

% save the cloud cover percentage
percent_cloud_cover = zeros(1, length(idx));


for nn = 1:length(idx)


    % -----------------------------------
    % ---- Write a Water Cloud file! ----
    % -----------------------------------
    % most uncertainties for the modis optical retrieval are between 2
    % and 10 percent. So lets round off all re values to the 1000th decimal
    % place

    % define the total cloud optical depth
    tau_c(nn) = round(modis.cloud.optThickness17(row(nn),col(nn)), 3);
    % define the droplet size

    re(nn) = round(modis.cloud.effRadius17(row(nn),col(nn)), 3);                      % microns

    % define the geometric location of the cloud top and cloud bottom
    if use_MODIS_cloudTopHeight==false
        z_topBottom = [9,8];          % km above surface

    else

        % if the cloud top height is below 1 km, make the lower altitude 0
        cloudTopHeight = modis.cloud.topHeight(row(nn), col(nn));

        if cloudTopHeight>=1000
            z_topBottom(nn,:) = [cloudTopHeight, cloudTopHeight - 1000]./1000;      % km above surface
        elseif cloudTopHeight<1000
            z_topBottom(nn,:) = [cloudTopHeight, 0]./1000;      % km above surface
        end

    end

    % ------------------------------------------------------
    % --------------------VERY IMPORTANT ------------------
    % ADD THE LOOP VARIABLE TO THE WC NAME TO MAKE IT UNIQUE
    % ------------------------------------------------------
    wc_filename{nn} = write_wc_file(re(nn),tau_c(nn),z_topBottom(nn,:), lambda_forTau, distribution_str,...
        distribution_var, vert_homogeneous_str, parameterization_str, nn);
    wc_filename{nn} = wc_filename{nn}{1};


    % ------------------------------------------------
    % ---- Define the input and output filenames! ----
    % ------------------------------------------------
    % input_names need a unique identifier. Let's give them the nn value so
    % they can be traced, and are writen over in memory
    if yesCloud==true
        inputName{nn} = [num2str(floor((wavelength(2)-wavelength(1))/2 + wavelength(1))),...
            'nm_withCloudLayer_',num2str(nn),'nn_',num2str(modis.solar.zenith(row(nn),col(nn))),'sza_',num2str(round(double(modis.sensor.zenith(row(nn),col(nn))))),...
            'vza_',atm_file(1:end-4),'.INP'];
    else
        inputName{nn} = [num2str((wavelength(2)-wavelength(1))/2 + wavelength(1)),...
            'nm_',num2str(nn),'nn_',num2str(modis.solar.zenith(row(nn),col(nn))),'sza_',num2str(round(double(modis.sensor.zenith(row(nn),col(nn))))),...
            'vza_',atm_file(1:end-4),'.INP'];
    end


    outputName{nn} = ['OUTPUT_',inputName{nn}(1:end-4)];



    % ----------------- ******************** ---------------------
    % ------------------ Write the INP File --------------------
    % ----------------- ******************** ---------------------

    % Open the old file for writing
    fileID = fopen([folder2save,inputName{nn}], 'w');

    % Define which RTE solver to use
    % ------------------------------------------------
    formatSpec = '%s %s %5s %s \n';
    fprintf(fileID, formatSpec,'rte_solver','disort',' ', '# Radiative transfer equation solver');


    % Define the number of streams to keep track of when solving the equation
    % of radiative transfer
    % ------------------------------------------------
    formatSpec = '%s %u %5s %s \n\n';
    fprintf(fileID, formatSpec,'number_of_streams', num_streams,' ', '# Number of streams');


    % Use phase function correction?
    % ------------------------------------------------
    if use_nakajima_phaseCorrection==true
        % define the pahse correction to be true
        % ------------------------------------------------
        formatSpec = '%s %5s %s \n\n';
        fprintf(fileID, formatSpec,'disort_intcor phase', ' ', '# Apply the Nakajima and Tanka radiance correction');
    end


    % Define the location and filename of the atmopsheric profile to use
    % ------------------------------------------------
    formatSpec = '%s %5s %s \n';
    fprintf(fileID, formatSpec,['atmosphere_file ','../data/atmmod/',atm_file],' ', '# Location of atmospheric profile');

    % Define the location and filename of the extraterrestrial solar source
    % ---------------------------------------------------------------------
    formatSpec = '%s %s %5s %s \n\n';
    fprintf(fileID, formatSpec,'source solar','../data/solar_flux/kurudz_1.0nm.dat', ' ', '# Bounds between 250 and 10000 nm');


    % Define the location and filename of the extraterrestrial solar source
    % ---------------------------------------------------------------------
    formatSpec = '%s %u %5s %s \n\n';
    fprintf(fileID, formatSpec,'day_of_year', day_of_year, ' ', '# accounts for changing Earth-Sun distance');


    % Define the total precipitable water
    % -------------------------------------------------------------------------
    % Define override of total precipitable water. This will force the total
    % column of water vapor to be whatever value you define.
    % If you don't wish to change the default, define the variable with nan
    total_h2O_column = modis.vapor.col_nir(row(nn),col(nn))*10;        % mm of precipitable water
    if isnan(total_h2O_column)==false
        formatSpec = '%s %s %f %s %5s %s \n';
        fprintf(fileID, formatSpec,'mol_modify','H2O', total_h2O_column,' MM', ' ', '# Total Precipitable Water');
    end


    % Define the surface albedo
    % ------------------------------------------------
    formatSpec = '%s %s %5s %s \n\n';
    fprintf(fileID, formatSpec,'albedo', albedo, ' ', '# Surface albedo of the ocean');

    
    % Define the Water Cloud properties, if you want a cloud in your model
    % --------------------------------------------------------------------
    if yesCloud==true

        % Define the water cloud file
        % ------------------------------------------------
        formatSpec = '%s %s %5s %s \n';
        %fprintf(fileID, formatSpec,'wc_file 1D', ['../data/wc/',wc_filename{nn}], ' ', '# Location of water cloud file');
        fprintf(fileID, formatSpec,'wc_file 1D', ['../data/wc/wc2/',wc_filename{nn}], ' ', '# Location of water cloud file');

        % Define the percentage of horizontal cloud cover
        % This is a number between 0 and 1
        % ------------------------------------------------
        if linear_cloudFraction==true
            percent_cloud_cover(nn) = -0.1082 * modis.EV1km.reflectance(row(nn), col(nn),1) + 0.92164;      % custom linear fit 1 (seventh pass figure)
            percent_cloud_cover(nn) = (0.89 - 0.8459)/(0.2 - 0.7) * modis.EV1km.reflectance(row(nn), col(nn),1) +...
                                        (0.8459 - 0.7*(0.89 - 0.8459)/(0.2 - 0.7));      % custom linear fit 2 (Ninth pass figure)
        else
            percent_cloud_cover(nn) = cloud_cover;
        end
        formatSpec = '%s %f %5s %s \n';
        fprintf(fileID, formatSpec,'cloudcover wc', percent_cloud_cover(nn), ' ', '# Cloud cover percentage');


        % Define the technique or parameterization used to convert liquid cloud
        % properties of r_eff and LWC to optical depth
        % ----------------------------------------------------------------------
        formatSpec = '%s %s %5s %s \n\n';
        fprintf(fileID, formatSpec,'wc_properties', wc_parameterization, ' ', '# optical properties parameterization technique');

    end

    % Define the wavelengths for which the equation of radiative transfer will
    % be solve
    % -------------------------------------------------------------------------
    formatSpec = '%s %f %f %5s %s \n\n';
    fprintf(fileID, formatSpec,'wavelength', wavelength(1), wavelength(2), ' ', '# Wavelength range');


    if use_coxMunk==true

        % Define the wind speed for the Cox-Munk ocean surface bi-directional reflectance model
        % be solve
        % -------------------------------------------------------------------------
        formatSpec = '%s %f %5s %s \n\n';
        fprintf(fileID, formatSpec,'brdf_cam u10', wind_speed, ' ', '# (m/s) Ocean Surface wind speed');

    end



    % Define the Aerosol Layer properties, if you want a cloud in your model
    % --------------------------------------------------------------------
    if yesAerosols==true

        % Turn on default aersol layer, which occupies lower 2km of model
        % --------------------------------------------------------------
        formatSpec = '%s %5s %s \n';
        fprintf(fileID, formatSpec,'aerosol_default', ' ', '# turn on Shettle (1989) boundary layer aerosols');


        % Specify the Aerosl type
        % 1=rural aersols,  4=maritime aersols,  5=Urban aerosols,
        % 6=Tropospheric aerosols
        % ------------------------------------------------
        formatSpec = '%s %u %5s %s \n';
        fprintf(fileID, formatSpec,'aerosol_haze', aerosol_type, ' ', '# Aerosol type');


        % Define aerosol layer optical depth
        % ----------------------------------------------------------------------
        formatSpec = '%s %f %5s %s \n\n';
        fprintf(fileID, formatSpec,'aerosol_modify tau set', aerosol_opticalDepth, ' ', '# Optical Depth of aerosol layer');

    end




    % Define the sensor altitude
    % ------------------------------------------------
    formatSpec = '%s %s %5s %s \n';
    fprintf(fileID, formatSpec,'zout', 'toa', ' ', '# Sensor Altitude');

    % Define the solar zenith angle
    % ------------------------------------------------
    % define the solar zenith angle
    sza = modis.solar.zenith(row(nn),col(nn));           % degree
    formatSpec = '%s %f %5s %s \n';
    fprintf(fileID, formatSpec,'sza', sza, ' ', '# Solar zenith angle');

    % Define the solar azimuth angle
    % -------------------------------------------------------
    % this is how we map MODIS azimuth of the sun to the LibRadTran measurement
    phi0 = modis.solar.azimuth(row(nn),col(nn)) + 180;         % degree
    formatSpec = '%s %f %5s %s \n';
    fprintf(fileID, formatSpec,'phi0', phi0, ' ', '# Solar azimuth angle');

    % Define the cosine of the zenith viewing angle
    % ------------------------------------------------
    % define the viewing zenith angle
    vza = double(modis.sensor.zenith(row(nn),col(nn))); % values are in degrees;                        % degree
    formatSpec = '%s %f %5s %s \n';
    fprintf(fileID, formatSpec,'umu', round(cosd(vza),4), ' ', '# Cosine of the zenith viewing angle');

    % Define the azimuth viewing angle
    % ------------------------------------------------
    % define the viewing azimuth angle
    % to properly map the azimuth angle onto the reference plane used by
    % libRadTran, we need an if statement
    if modis.sensor.azimuth(row(nn),col(nn))<0
        vaz = 360 + modis.sensor.azimuth(row(nn),col(nn));
    else
        vaz = modis.sensor.azimuth(row(nn),col(nn));
    end
    formatSpec = '%s %f %5s %s \n\n';
    fprintf(fileID, formatSpec,'phi', vaz, ' ', '# Azimuthal viewing angle');


    % Set the error message to quiet of verbose
    % ------------------------------------------------
    formatSpec = '%s';
    fprintf(fileID, formatSpec,'verbose');


    % Close the file!
    fclose(fileID);
    % ----------------------------------------------------
    % ----------------------------------------------------


end



%% Run INP files

R_model = zeros(1, length(idx));
R_modis = zeros(1, length(idx));
uncertainty = zeros(1, length(idx));

% store the modis reflectance and uncertainties in a standalone array to
% avoid unecessary communication during the parallel for loop

modis_reflectance = modis.EV1km.reflectance(:, :, 1);
modis_reflectance_uncertainty = modis.EV1km.reflectanceUncert(:,:,1);

tic
parfor nn = 1:length(idx)

    % Store the modis reflectance value
     R_modis(nn) =  modis_reflectance(row(nn), col(nn));

     % store the MODIS reflectance uncertainty
     uncertainty(nn) = 0.01*modis_reflectance_uncertainty(row(nn), col(nn)) *  R_modis(nn);

    % compute INP file
    [inputSettings] = runUVSPEC(folder2save,inputName{nn},outputName{nn});

    % read .OUT file
    [ds,~,~] = readUVSPEC(folder2save,outputName{nn},inputSettings(2,:));

    % save reflectance
    R_model(nn) = reflectanceFunction(inputSettings(2,:), ds, ones(length(ds.wavelength),1));

end
toc

%% Plot the ratio of modeled reflectance with MODIS reflectance to the optical depth


figure; 
subplot(1,3,1)
plot(tau_c, R_model./R_modis, '.', 'MarkerSize', 25)


hold on; grid on; grid minor

xlabel('$\tau_c$','Interpreter','latex', 'FontSize',35)
ylabel('$\frac{R_{est}}{R_{modis}}$','Interpreter','latex', 'FontSize',35)
title('Ratio doesnt trend well with optical depth', 'FontSize',25)
set(gcf, 'Position', [0 0 2000 800])
axis square

subplot(1,3,2)
plot(R_modis, R_model./R_modis, '.', 'MarkerSize', 25)


hold on; grid on; grid minor

xlabel('Measured Reflectance $(1/sr)$','Interpreter','latex', 'FontSize',35)
ylabel('$\frac{R_{est}}{R_{modis}}$','Interpreter','latex', 'FontSize',35)
title('Ratio doesnt trend well with measured reflectance', 'FontSize',25)
axis square

subplot(1,3,3)
plot(R_model, R_modis, '.', 'MarkerSize', 25)


hold on; grid on; grid minor

xlabel('My Estimate of Reflectance $(1/sr)$','Interpreter','latex', 'FontSize',25)
ylabel('MODIS Measured Reflectance $(1/sr)$','Interpreter','latex', 'FontSize',25)
title('Deviation from 1 to 1', 'FontSize',25)
axis square




%% Sub Plot comparing modeled and measured reflectance, re, tau_c, and a histogram

label_fontSize = 20;
% ------------------------ SUBPLOT 1 ----------------------------------
% First make a scatter plot comparing the modeled and measured
% reflectance and have the color of each marker represent the retrieved
% droplet size
% ---------------------------------------------------------------------


% create 1 to 1 line
% find the minimum and maximum values to create a y=x line

min_R_est = min(R_model);
min_R_modis = min(R_modis);

max_R_est = max(R_model);
max_R_modis = max(R_modis);

min_R_global = min([min_R_est,min_R_modis]);

max_R_global = max([max_R_est,max_R_modis]);

x_r = linspace((0.9 * min_R_global),(1.1*max_R_global),150);

% Lets define the color of each marker to be associated with the droplet
% size
% set the number of colors to be the length of the data to plot
C = colormap(parula(length(R_model)));
% sort the droplet size values
[re_sort, idx_sort] = sort(re, 'ascend');

figure; subplot(1,3,1)
plot(x_r, x_r, 'k', 'LineWidth',1)
hold on

for nn = 1:length(re_sort)
    
   errorbar(R_model(idx_sort(nn)), R_modis(idx_sort(nn)), uncertainty(idx_sort(nn)),'vertical','Marker','.','Color',C(nn,:),'MarkerSize',25)

   hold on

end

% set the colorbar limits
% set the limits of the colormap to be the min and max value
cb = colorbar;
clim([min(re), max(re)]);
% set colorbar title
cb.Label.String = '$r_e$ ($\mu m$)';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = 25;

hold on; grid on; grid minor
xlim([x_r(1), x_r(end)])
ylim([x_r(1), x_r(end)])
xlabel('My Estimate of Reflectance $(1/sr)$','Interpreter','latex', 'FontSize',label_fontSize)
ylabel('MODIS Measured Reflectance $(1/sr)$','Interpreter','latex', 'FontSize',label_fontSize)
set(gcf, 'Position', [0 0 1800 900])
axis square


% ------------------------ SUBPLOT 2 ----------------------------------
% Next make a scatter plot comparing the modeled and measured
% reflectance and have the color of each marker represent the retrieved
% optical depth
% ---------------------------------------------------------------------



% Lets define the color of each marker to be associated with the droplet
% size
% set the number of colors to be the length of the data to plot
C = colormap(parula(length(tau_c)));
% sort the droplet size values
[tau_c_sort, idx_sort] = sort(tau_c, 'ascend');

subplot(1,3,2)
plot(x_r, x_r, 'k', 'LineWidth',1)
hold on

for nn = 1:length(tau_c_sort)
    
   errorbar(R_model(idx_sort(nn)), R_modis(idx_sort(nn)), uncertainty(idx_sort(nn)),'vertical','Marker','.','Color',C(nn,:),'MarkerSize',25)

   hold on

end

% set the colorbar limits
% set the limits of the colormap to be the min and max value
cb = colorbar;
clim([min(tau_c), max(tau_c)]);
% set colorbar title
cb.Label.String = '$\tau_c$ ($\mu m$)';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = 25;

hold on; grid on; grid minor
xlim([x_r(1), x_r(end)])
ylim([x_r(1), x_r(end)])
xlabel('My Estimate of Reflectance $(1/sr)$','Interpreter','latex', 'FontSize',label_fontSize)
ylabel('MODIS Measured Reflectance $(1/sr)$','Interpreter','latex', 'FontSize',label_fontSize)
axis square



% ------------------------ SUBPLOT 3 ----------------------------------
% Make a histogram of the ratio between my reflectance estiamte and the
% true measured reflectance
% ---------------------------------------------------------------------


subplot(1,3,3)

h = histogram(R_model./R_modis, 'NumBins',50);
% find the bin with the most counts - the mode
[max_count, idx_max] = max(h.BinCounts);
mode_val = h.BinEdges(idx_max) + h.BinWidth/2;
% let's also compute the weighted average
avg_weighted = sum((h.BinCounts./sum(h.BinCounts)) .* (h.BinEdges(1:end-1) + h.BinWidth/2))...
                / sum((h.BinCounts./sum(h.BinCounts)));


hold on
% Plot the modal value
% xline(mode_val,'Label',['Mode = ',num2str(round(mode_val,2))], 'FontSize', 15,...
%     'LineWidth',1)
% Plot the weighted average
xline(mode_val,'Label',['W-Mean = ',num2str(round(avg_weighted,2))], 'FontSize', 17,...
    'LineWidth',1.5, 'FontWeight','bold', 'LabelVerticalAlignment','bottom')

hold on; grid on; grid minor

xlabel('$\frac{R_{est}}{R_{modis}}$','Interpreter','latex', 'FontSize',30)
ylabel('Counts','Interpreter','latex', 'FontSize',label_fontSize)
axis square
