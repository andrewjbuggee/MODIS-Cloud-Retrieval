%% ----- READ IN L2 MODIS DATA -----

% this function will read in L2 MODIS data as .hdf files
% will produce all necessary information into a cell array, and the data
% set into a structure

% there are many different fields of data one could read from a MODIS .hdf


% ---- Description of the different fields within the HDF File -----




% By Andrew J. Buggee
%%

function [cloud] = readMODIS_L2_data(fileName)
%% ---- Read in Conversion Scales and Offsets -----

% retreive the scales and offsets for effective particle radius and optical
% thickness

cloudProp_info = hdfinfo(fileName);


% -----------------------------------------------------
% -------- PULL SCALE FACTORS AND OFFSETS -------------
% -----------------------------------------------------

% extract effective radius info first

effectiveRadius17_scales = cloudProp_info.Vgroup.Vgroup(2).SDS(67).Attributes(5).Value; %  
effectiveRadius17_offset = cloudProp_info.Vgroup.Vgroup(2).SDS(67).Attributes(6).Value;

% effective radius uncertainty for bands 1 and 7
effectiveRadius_uncertainty_17_scales = cloudProp_info.Vgroup.Vgroup(2).SDS(91).Attributes(5).Value;
effectiveRadius_uncertainty_17_offset = cloudProp_info.Vgroup.Vgroup(2).SDS(91).Attributes(6).Value;

% Effective radius scale factor and offset for bands 1 and 6
effectiveRadius16_scales = cloudProp_info.Vgroup.Vgroup(2).SDS(69).Attributes(5).Value; %  
effectiveRadius16_offset = cloudProp_info.Vgroup.Vgroup(2).SDS(69).Attributes(6).Value;


% effective radius uncertainty for bands 1 and 6
% uncertainties are listed as percents
effectiveRadius_uncertainty_16_scales = cloudProp_info.Vgroup.Vgroup(2).SDS(92).Attributes(5).Value;
effectiveRadius_uncertainty_16_offset = cloudProp_info.Vgroup.Vgroup(2).SDS(92).Attributes(6).Value;



% extract the optical thickness scales and offset
optThickness17_scales = cloudProp_info.Vgroup.Vgroup(2).SDS(73).Attributes(5).Value; % 
optThickness17_offset = cloudProp_info.Vgroup.Vgroup(2).SDS(73).Attributes(6).Value;


% optical thicknee uncertainty for bands 1 and 7
% uncertainties are listed as percents
optThickness_uncertainty_17_scales = cloudProp_info.Vgroup.Vgroup(2).SDS(94).Attributes(5).Value;
optThickness_uncertainty_17_offset = cloudProp_info.Vgroup.Vgroup(2).SDS(94).Attributes(6).Value;


% optical thickness scale factors and offsets for bands 1 and 6
optThickness16_scales = cloudProp_info.Vgroup.Vgroup(2).SDS(75).Attributes(5).Value; %  output will be a cell array
optThickness16_offset = cloudProp_info.Vgroup.Vgroup(2).SDS(75).Attributes(6).Value;

% optical thicknee uncertainty for bands 1 and 6
% uncertainties are listed as percents
optThickness_uncertainty_16_scales = cloudProp_info.Vgroup.Vgroup(2).SDS(95).Attributes(5).Value;
optThickness_uncertainty_16_offset = cloudProp_info.Vgroup.Vgroup(2).SDS(95).Attributes(6).Value;

% extract the cloud top height at 1km resolution scales and offset
cloudTopHeight_scales = cloudProp_info.Vgroup.Vgroup(2).SDS(58).Attributes(5).Value; %  output will be a cell array
cloudTopHeight_offset = cloudProp_info.Vgroup.Vgroup(2).SDS(58).Attributes(6).Value;

% extract the cloud top pressure at 1km resolution scales and offset
cloudTopPressure_scales = cloudProp_info.Vgroup.Vgroup(2).SDS(57).Attributes(5).Value; %  output will be a cell array
cloudTopPressure_offset = cloudProp_info.Vgroup.Vgroup(2).SDS(57).Attributes(6).Value;

% extract the cloud top temperature at 1km resolution scales and offset
cloudTopTemperature_scales = cloudProp_info.Vgroup.Vgroup(2).SDS(59).Attributes(5).Value; %  output will be a cell array
cloudTopTemperature_offset = cloudProp_info.Vgroup.Vgroup(2).SDS(59).Attributes(6).Value;

% extract the cloud phase used in Optical Thickness/Effective Radius determination -  scales and offset
% The values in this SDS are set to mean the following:                              
% 0 -- cloud mask undetermined                                                       
% 1 -- clear sky                                                                     
% 2 -- liquid water cloud                                                            
% 3 -- ice cloud                                                                     
% 4 -- undetermined phase cloud (but retrieval is attempted as  liquid water)        

cloudPhase_scales = cloudProp_info.Vgroup.Vgroup(2).SDS(105).Attributes(5).Value; %  output will be a cell array
cloudPhase_offset = cloudProp_info.Vgroup.Vgroup(2).SDS(105).Attributes(6).Value;




% -----------------------------------------------------
% ---------------- PULL THE DATA SET ------------------
% -----------------------------------------------------inties are listed as percents


% extract the effective radius data using bands 1 and 7
effectiveRadius_17 = hdfread(fileName,'Cloud_Effective_Radius');

% extract the effective radius uncertainty bands 1 and 7
% uncertainties are listed as percents
effectiveRadius_uncertainty_17 = hdfread(fileName,'Cloud_Effective_Radius_Uncertainty');

% extract the effective radius data using bands 1 and 6
effectiveRadius_16 = hdfread(fileName,'Cloud_Effective_Radius_16');

% extract the effective radius uncertainty bands 1 and 6
% uncertainties are listed as percents
effectiveRadius_uncertainty_16 = hdfread(fileName,'Cloud_Effective_Radius_Uncertainty_16');


% extract the Optical thickness data using bands 1 and 7
opticalThickness_17 = hdfread(fileName,'Cloud_Optical_Thickness');

% extract the optical thickness uncertainty bands 1 and 7
% uncertainties are listed as percents
optThickness_uncertainty_17 = hdfread(fileName,'Cloud_Optical_Thickness_Uncertainty');



% extract the Optical thickness data using bands 1 and 6
opticalThickness_16 = hdfread(fileName,'Cloud_Optical_Thickness_16');

% extract the optical thickness uncertainty bands 1 and 7
% uncertainties are listed as percents
optThickness_uncertainty_16 = hdfread(fileName,'Cloud_Optical_Thickness_Uncertainty_16');


% extract the cloud top geopotential height at 1km resolution - units (meters)
cloudTop_geopotentialHeight = hdfread(fileName,'cloud_top_height_1km');

% extract the cloud top pressure - units (hpa)
cloudTop_pressure = hdfread(fileName,'cloud_top_pressure_1km');

% extract the cloud top temperature - units (k)
cloudTop_temperature = hdfread(fileName,'cloud_top_temperature_1km');

% extract the cloud phase - 
% The values are:
%   0 - no phase result
%   1 - no phase result
%   2 - liquid water
%   3 - ice
%   4 - undetermined

cloudPhase = hdfread(fileName,'Cloud_Phase_Optical_Properties');


% Read the cloud_mask_SPI variable which contains the sub-pixel
% heterogeneity index. Larger values have been shown the have more
% retrieval bias and increased rates of retreival failure

subPix_heteroIndex = hdfread(fileName, 'Cloud_Mask_SPI');

subPix_heteroIndex_scales = cloudProp_info.Vgroup.Vgroup(2).SDS(118).Attributes(5).Value;
subPix_heterIndex_offset = cloudProp_info.Vgroup.Vgroup(2).SDS(118).Attributes(6).Value;

%% --- Convert Data ---


cloud.effRadius17 = scalesOffsets2Matrix(effectiveRadius_17,effectiveRadius17_scales,effectiveRadius17_offset);
cloud.effRad_uncert_17 = scalesOffsets2Matrix(effectiveRadius_uncertainty_17,effectiveRadius_uncertainty_17_scales,effectiveRadius_uncertainty_17_offset);
cloud.optThickness17 = scalesOffsets2Matrix(opticalThickness_17,optThickness17_scales,optThickness17_offset);
cloud.optThickness_uncert_17 = scalesOffsets2Matrix(optThickness_uncertainty_17,optThickness_uncertainty_17_scales,optThickness_uncertainty_17_offset);

cloud.effRadius16 = scalesOffsets2Matrix(effectiveRadius_16,effectiveRadius16_scales,effectiveRadius16_offset);
cloud.effRad_uncert_16 = scalesOffsets2Matrix(effectiveRadius_uncertainty_16,effectiveRadius_uncertainty_16_scales,effectiveRadius_uncertainty_16_offset);
cloud.optThickness16 = scalesOffsets2Matrix(opticalThickness_16,optThickness16_scales,optThickness16_offset);
cloud.optThickness_uncert_16 = scalesOffsets2Matrix(optThickness_uncertainty_16,optThickness_uncertainty_16_scales,optThickness_uncertainty_16_offset);

cloud.topHeight = scalesOffsets2Matrix(cloudTop_geopotentialHeight,cloudTopHeight_scales,cloudTopHeight_offset);
cloud.topPressure = scalesOffsets2Matrix(cloudTop_pressure,cloudTopPressure_scales,cloudTopPressure_offset);
cloud.topTemperature = scalesOffsets2Matrix(cloudTop_temperature,cloudTopTemperature_scales,cloudTopTemperature_offset);

cloud.phase = scalesOffsets2Matrix(cloudPhase,cloudPhase_scales,cloudPhase_offset);

cloud.SPI = scalesOffsets2Matrix(subPix_heteroIndex,subPix_heteroIndex_scales,subPix_heterIndex_offset);







end





