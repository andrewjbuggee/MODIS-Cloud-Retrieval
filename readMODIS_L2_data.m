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

% extract effective radius info first

effectiveRadius17_scales = cloudProp_info.Vgroup.Vgroup(2).SDS(67).Attributes(5).Value; %  output will be a cell array
effectiveRadius17_offset = cloudProp_info.Vgroup.Vgroup(2).SDS(67).Attributes(6).Value;

% effective radius uncertainty for bands 1 and 7
effectiveRadius_uncertainty_17_scales = cloudProp_info.Vgroup.Vgroup(2).SDS(91).Attributes(5).Value;
effectiveRadius_uncertainty_17_offset = cloudProp_info.Vgroup.Vgroup(2).SDS(91).Attributes(6).Value;

% effective radius uncertainty for bands 1 and 6
effectiveRadius_uncertainty_16_scales = cloudProp_info.Vgroup.Vgroup(2).SDS(92).Attributes(5).Value;
effectiveRadius_uncertainty_16_offset = cloudProp_info.Vgroup.Vgroup(2).SDS(92).Attributes(6).Value;


effectiveRadius16_scales = cloudProp_info.Vgroup.Vgroup(2).SDS(69).Attributes(5).Value; %  output will be a cell array
effectiveRadius16_offset = cloudProp_info.Vgroup.Vgroup(2).SDS(69).Attributes(6).Value;

% extract the optical thickness scales and offset
optThickness17_scales = cloudProp_info.Vgroup.Vgroup(2).SDS(73).Attributes(5).Value; %  output will be a cell array
optThickness17_offset = cloudProp_info.Vgroup.Vgroup(2).SDS(73).Attributes(6).Value;

optThickness16_scales = cloudProp_info.Vgroup.Vgroup(2).SDS(75).Attributes(5).Value; %  output will be a cell array
optThickness16_offset = cloudProp_info.Vgroup.Vgroup(2).SDS(75).Attributes(6).Value;

% extract the cloud top height scales and offset
cloudTopHeight_scales = cloudProp_info.Vgroup.Vgroup(2).SDS(18).Attributes(5).Value; %  output will be a cell array
cloudTopHeight_offset = cloudProp_info.Vgroup.Vgroup(2).SDS(18).Attributes(6).Value;

% extract the cloud top pressure scales and offset
cloudTopPressure_scales = cloudProp_info.Vgroup.Vgroup(2).SDS(22).Attributes(5).Value; %  output will be a cell array
cloudTopPressure_offset = cloudProp_info.Vgroup.Vgroup(2).SDS(22).Attributes(6).Value;

% extract the cloud top temperature scales and offset
cloudTopTemperature_scales = cloudProp_info.Vgroup.Vgroup(2).SDS(28).Attributes(5).Value; %  output will be a cell array
cloudTopTemperature_offset = cloudProp_info.Vgroup.Vgroup(2).SDS(28).Attributes(6).Value;

% extract the cloud top temperature scales and offset
cloudPhase_scales = cloudProp_info.Vgroup.Vgroup(2).SDS(51).Attributes(5).Value; %  output will be a cell array
cloudPhase_offset = cloudProp_info.Vgroup.Vgroup(2).SDS(51).Attributes(6).Value;


% extract the effective radius data using bands 1 and 7
effectiveRadius_17 = hdfread(fileName,'Cloud_Effective_Radius');

% extract the effective radius uncertainty bands 1 and 7
effectiveRadius_uncertainty_17 = hdfread(fileName,'Cloud_Effective_Radius_Uncertainty');

% extract the effective radius data using bands 1 and 6
effectiveRadius_16 = hdfread(fileName,'Cloud_Effective_Radius_16');

% extract the effective radius uncertainty bands 1 and 6
effectiveRadius_uncertainty_16 = hdfread(fileName,'Cloud_Effective_Radius_Uncertainty_16');

% extract the Optical thickness data using bands 1 and 7
opticalThickness_17 = hdfread(fileName,'Cloud_Optical_Thickness');

% extract the Optical thickness data using bands 1 and 6
opticalThickness_16 = hdfread(fileName,'Cloud_Optical_Thickness_16');

% extract the cloud top geopotential height - units (meters)
cloudTop_geopotentialHeight = hdfread(fileName,'Cloud_Top_Height');

% extract the cloud top pressure - units (hpa)
cloudTop_pressure = hdfread(fileName,'Cloud_Top_Pressure');

% extract the cloud top temperature - units (k)
cloudTop_temperature = hdfread(fileName,'Cloud_Top_Temperature');

% extract the cloud phase - 
% The values are:
%   0 - no phase result
%   1 - no phase result
%   2 - liquid water
%   3 - ice
%   4 - undetermined

cloudPhase = hdfread(fileName,'Cloud_Phase_Optical_Properties');

%% --- Convert Data ---


cloud.effRadius17 = scalesOffsets2Matrix(effectiveRadius_17,effectiveRadius17_scales,effectiveRadius17_offset);
cloud.effRad_uncertainty_17 = scalesOffsets2Matrix(effectiveRadius_uncertainty_17,effectiveRadius_uncertainty_17_scales,effectiveRadius_uncertainty_17_offset);
cloud.optThickness17 = scalesOffsets2Matrix(opticalThickness_17,optThickness17_scales,optThickness17_offset);

cloud.effRadius16 = scalesOffsets2Matrix(effectiveRadius_16,effectiveRadius16_scales,effectiveRadius16_offset);
cloud.effRad_uncertainty_16 = scalesOffsets2Matrix(effectiveRadius_uncertainty_16,effectiveRadius_uncertainty_16_scales,effectiveRadius_uncertainty_16_offset);
cloud.optThickness16 = scalesOffsets2Matrix(opticalThickness_16,optThickness16_scales,optThickness16_offset);

cloud.topHeight = scalesOffsets2Matrix(cloudTop_geopotentialHeight,cloudTopHeight_scales,cloudTopHeight_offset);
cloud.topPressure = scalesOffsets2Matrix(cloudTop_pressure,cloudTopPressure_scales,cloudTopPressure_offset);
cloud.topTemperature = scalesOffsets2Matrix(cloudTop_temperature,cloudTopTemperature_scales,cloudTopTemperature_offset);

cloud.phase = scalesOffsets2Matrix(cloudPhase,cloudPhase_scales,cloudPhase_offset);







end





