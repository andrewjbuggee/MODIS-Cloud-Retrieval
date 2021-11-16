%% ----- Estimate MODIS Cloud Optical Properties -----


clear variables;
% By Andrew J. Buggee

%% ----- Extract MODIS data of Interest -----


% Define data folders and files that you'd like to read

% add the libradtran directory to the path

% this one is for my LASP computer
% modisINP_folderName = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/MODIS_08_25_2021/';

% this one is for my personal laptop
modisINP_folderName = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/',...
    'Hyperspectral-Cloud-Droplet-Retrieval-Research/LibRadTran/libRadtran-2.0.4/MODIS_08_25_2021/'];

% define the files names

L1B_fileName = 'MOD02HKM.A2021237.1920.006.2021238072711.hdf';

geoFileName = 'MOD03.A2021237.1920.006.2021238012056.hdf';

L2_fileName = 'MOD06_L2.A2021237.1920.061.2021238075516.hdf';

modis = retrieveMODIS_data(L1B_fileName,geoFileName,L2_fileName);

%% ----- Create a structure defining inputs of the problem -----

% this is a built-in function that is defined at the bottom of this script
inputs = create_modis_inputs(L1B_fileName);

disp('Check inputs to make sure they are what you want!!')

%% ----- Find suitable Pixels! -----

% find pixels within the modis data set that fit our needs
% if we've alraedy done this long calculation, we just load the pixel set
% from a saved .mat file

% lets check to see if there is a suitable pixels .mat file in our folder

pixels_file_flag = isfile([inputs.savedCalculations_folderName,'suitablePixels.mat']);

if inputs.flags.findSuitablePixels == true && pixels_file_flag == false
    
    pixels = findSuitablePixel(modis,inputs);
    
    % save these pixels in a .mat file, along with the inputs
    save([inputs.savedCalculations_folderName,'suitablePixels.mat'],'pixels','inputs');
    
    % now we want only a subset of pixels to write INP files
    
    pixels2use = load_suitablePixels(inputs);
    
elseif pixels_file_flag == true && inputs.flags.loadPixelSet == false
    
    % we don't need to load the entire pixels file into out workspace. But
    % if we do load a subset of pixels, we need a way to trace back to what
    % pixels we used. Save the pixels used to the data folder listed in the
    % inputs
    
    pixels2use = load_suitablePixels(inputs);

elseif pixels_file_flag == true && inputs.flags.loadPixelSet == true
    
    disp('Load a pixel set from a previous calculation!')
    
end




%% ----- Create .INP files for MODIS Geometry -----


if inputs.flags.writeINPfiles == true
    % which pixels on the MODIS array are we using for the gemoetry of the
    % problem? The pixels that we found to be suitable!
    
    names.inp = write_INP_4_MODIS_hdf(inputs,pixels2use,modis);
    
    % now lets write the output names
    
    names.out = writeOutputNames(names.inp);
else
    
    % if the files already exist, just grab the names!
    names.inp = getMODIS_INPnames_withClouds(modis.solar,inputs);
    names.out = writeOutputNames(names.inp);
end

%% ----- Run uvspec and calculate Reflectance Function for Model -----

% geometry stays the same, but we calculate the radiative transfer equation
% for different cloud values (tau,re)

if inputs.flags.runUVSPEC == true
    
    % R is the reflectance integrated over a bandwidth
    % Rl is the reflectance at each spectral bin
    [R,Rl] = runReflectanceFunction(inputs,names);
    
elseif inputs.flags.runUVSPEC == false
    
    
    
end
%% ----- Compute the Reflectance Function for the MODIS Observations -----

% We don't have to calculate the reflectance function of MODIS if we don't
% want to. They proivde it as an output in their data

modisR = grab_modis_reflectance(modis,inputs);


%% ----- Compare Reflectance Fucntion of MODIS with Theoretical Calculations (Grid Search) -----

% first grid search is on a coarse grid
% we want to minimize two the reflectance for two wavelengths

% if interpGridScalFactor is 10, then 9 rows will be interpolated to be 90
% rows, and 10 columns will be interpolated to be 100 columns

minVals = leastSquaresGridSearch(modisR, R, inputs, pixels2use);

[teTable, reflTable] = gatherTruthEstimateVals(modis, minVals, inputs); % containts truth ad estimates and difference


%% ----- Makes Plots -----

% Plot MODIS Calculated cloud properties
figure; subplot(1,2,1); imagesc(modis.cloud.effRadius); colorbar; title('effective radius (\mum)')
subplot(1,2,2); imagesc(modis.cloud.optThickness); colorbar; title('optical thickness')


% plot model reflectance function

plot2ReflectanceFuncBands(modis,R,inputs)

%%
















%% ----- CREATE INPUTS -----

function inputs = create_modis_inputs(L1B_fileName)


% ----- what water cloud files should we use? -----

% the values for re and tau will be the values that are used in the inp
% files for the liquid water cloud cloud files. These are already written.
% Should probably fix this later

inputs.re = [4:4:20,25:5:40]; % - microns - effective radius - value is in the file name
inputs.tau_c = [1,5:5:30,40:10:80]; % cloud optical thickness - value is in the file name



inputs.bands2run = [1,2,6,7]; % these are the bands that we will run uvspec with
inputs.bands2search = [1,7]; % these are the modis bands that are used in the retrieval problem
inputs.bands2plot = [1,7]; % these are the modis bands that will be plotted, both the modis calcualted stuff and the stuff I calcualte

% if interpGridScalFactor is 10, then 9 rows will be interpolated to be 90
% rows, and 10 columns will be interpolated to be 100 columns
inputs.interpGridScaleFactor = 10; % scale factor the will be used to increase the grid size for interpolation.

inputs.savedCalculations_folderName = './MODIS_data/2021_08_25/'; % this is the folder that all the saved calculations for the 2021_08_25 data set will go
inputs.saveCalculations_fileName = ['uvspec_CALCS_',date,'.mat'];

inputs.INP_folderName = ['MODIS_day_',L1B_fileName(15:17),'_year_',L1B_fileName(11:14),...
    '_time_',L1B_fileName(19:22),'/']; % this is the folder name that the INP files will be written to 


% ----- defining what pixels to use -----

inputs.pixels.tauThreshold = 15; % only find pixels in the modis data that is greater than or equal to a tau of this 
inputs.pixels.num_INP_files_2write = 1; % we will randomly select this many pixels from the set of suitable pixels found to create .INP files




% ----- flags -----

% define flags that tell the codes to either run certain things, or don't
% run certain things

inputs.flags.findSuitablePixels = false; % if true, this will search the modis data set for pixels to use
inputs.flags.loadPixelSet = true; % if true, the code will load an older set of pixels that has already been used before, and likely has INP files. If false, it tells the code to find a new subset of pixels
inputs.flags.writeINPfiles = false; % if true, this will create inp files for each the length of vector pixel.row
inputs.flags.runUVSPEC = false; % if true, this will run all of the inp files create from the above flag through uvspec
inputs.flags.plotMLS_figures = false; % this will tell the leasSquaresGridSearch code to plot the l



end

