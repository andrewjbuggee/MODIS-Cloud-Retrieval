%% ----- CREATE INPUTS NEEDED TO COMPUTE TBLUT METHOD ON MODIS DATA -----


% INPUTS:
%   (1) folderName - 

%   (2) L1B_fileName - 


% OUTPUTS:
%   (1) inputs - effective droplet radius profile


% By Andrew John Buggee
%%

function inputs = create_modis_inputs(folderName, L1B_fileNames)


% --- SAVE THE MODIS FILE NAME ----
inputs.modisDataFolder = folderName;

% ----- what water cloud files should we use? -----

% The retrieved values of re and Tau computed by the MODIS TBLUT algorithm
% will be used to write liquid water cloud files for LibRadTran. 


% ?? These are already written. Should probably fix this later ??

% These are the values that the forward model will use to create a
% pre-computed table of reflectances for the specfic geometry of each pixel

inputs.re = [1:3:25]; % - microns - effective radius - value is in the file name
inputs.tau_c = [1,5:5:80]; % cloud optical thickness - value is in the file name



inputs.bands2run = [1,6,7]; % these are the bands that we will run uvspec with
inputs.bands2search = [1,7; 1,6]; % these are the modis bands that are used in the retrieval problem
inputs.bands2plot = [1,7]; % these are the modis bands that will be plotted, both the modis calcualted stuff and the stuff I calcualte

% if interpGridScaleFactor is 10, then 9 rows will be interpolated to be 90
% rows, and 10 columns will be interpolated to be 100 columns
inputs.interpGridScaleFactor = 150; % scale factor the will be used to increase the grid size for interpolation.

inputs.savedCalculations_folderName = folderName; % this is the folder that all the saved calculations will go
inputs.saveCalculations_fileName = ['uvspec_CALCS_',date,'.mat'];

inputs.INP_folderName = ['MODIS_day_',L1B_fileNames{1}(15:17),'_year_',L1B_fileNames{1}(11:14),...
    '_time_',L1B_fileNames{1}(19:22),'/']; % this is the folder name that the INP files will be written to 


% ----- defining what pixels to use -----

% only find pixels in the modis data that is greater than or equal to a tau 
% defined by the value below
inputs.pixels.tauThreshold = 2; 
 
% we will randomly select this many pixels from the set of suitable pixels
% found to create .INP files Each pixel and its associated geometry will
% have an INP file of its own that will solve the equation of radiative
% transfer
inputs.pixels.num_2calculate = 100;


% ------------------
% ----- FLAGS! -----
% ------------------

% define flags that tell the codes to either run certain things, or don't
% run certain things

inputs.flags.findSuitablePixels = false; % if true, this will search the modis data set for pixels to use

% if true, the code will load an older set of pixels that has already been used before, and 
% likely has INP files. If false, it tells the code to find a new random subset of pixels
inputs.flags.loadPixelSet = false; 
inputs.flags.writeINPfiles = true; % if true, this will create inp files for each the length of vector pixel.row
inputs.flags.runUVSPEC = true; % if true, this will run all of the inp files create from the above flag through uvspec
inputs.flags.plotMLS_figures = false; % this will tell the leasSquaresGridSearch code to plot the 


% ------------------
% ----- Stuff for writing water cloud files! -----
% ------------------
% can be 'hu' or 'mie interpolate'
inputs.flags.wc_properties = 'mie interpolate';        % use the hu and stamnes parameterization for converting cloud properties to optical properties
% can either be 'mie' or '2limit'
inputs.flags.wc_parameterization = 'mie';        % This string is used to compute the LWC from optical depth and effective radius
% can either be 'mono' or 'gamma'
inputs.flags.distribution_type = 'gamma';        % This string is used to compute the LWC from optical depth and effective radius

% ----- ISSUE A WARNING! SETTINGS SHOULD BE CHECKED -----

warning([newline, 'Check inputs structure to make sure the settings reflect the situation you wish to model!', newline]);

end
