%% Plot the VOCALS-REX flight path and the overlapping MODIS pixels

% This will plot all of the MODIS pixels found using modisIndex_min_dist
% These are the unique pixels found closest to the VOCALS-REx profile.


% By Andrew John Buggee

%%

function plot_MODIS_pixels_overlapping_VocalsRex_profile(modis,vocalsRex, modisFolder)

 

figure; 


% load the across and along pixel growth curve fits
if strcmp(whatComputer, 'anbu8374')==true

    load(['/Users/anbu8374/Documents/MATLAB/HyperSpectral_Cloud_Retrieval/MODIS_Cloud_Retrieval/',...
        'along_across_pixel_growth.mat'])

elseif strcmp(whatComputer, 'andrewbuggee')==true

    load(['/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/',...
        'MODIS_Cloud_Retrieval/along_across_pixel_growth.mat'])

end

% Assume these are the center coordinates of the pixel
% First, create a perfect square assuming, when MODIS is looking nadir,
% that the pixel is 1km by 1 km on the ground.
% Thus, the diagonal from the center to any corner is 1/sqrt(2) km long
%dist_to_corners = 1/sqrt(2);        % km

% Using the 'reckon' function, we con move to each corner from the center
% position. But we need to know the orientation of the MODIS instrument
% with respect to the meriodional lines.

terra_inclination = 98.2098;        % degrees
aqua_inclination = 	98.1987;        % degrees

azimuth_2_corners = [45, 135, 225, 315, 45];


% grab the indexes of the ordered set of effective radius for the data
% being used
data2plot = modis.cloud.effRadius17(vocalsRex.modisIndex_minDist);
% Set NaNs to 0
data2plot(isnan(data2plot))=0;

[~, index_sort] = sort(reshape(data2plot, [],1), 'ascend');
%[r,c] = ind2sub([n_rows, n_cols], index_sort);
C = parula(length(index_sort));


for ii = 1:length(data2plot)

        % the distance to the corners depends on the viewing zenith angle
        along_length = along_scan(modis.sensor.zenith(vocalsRex.modisIndex_minDist(ii)));
        across_length = across_scan(modis.sensor.zenith(vocalsRex.modisIndex_minDist(ii)));

        % compute the distance from the center to each corner
        dist_to_corner = sqrt((along_length/2)^2 + (across_length/2)^2) * 1e3;  % meters


        % The azimuth angle in the reckon function is with respect to the local
        % meridian (north). So we have to add the inclination angle to the azimuth
        % angles listed above.

        [lat_corner, long_corner] = reckon(modis.geo.lat(vocalsRex.modisIndex_minDist(ii)),...
            modis.geo.long(vocalsRex.modisIndex_minDist(ii)),linspace(dist_to_corner, dist_to_corner,5),...
            azimuth_2_corners+terra_inclination, wgs84Ellipsoid);

        % Create a geopolyshape
        modis_polyshape = geopolyshape(lat_corner, long_corner);

        gp = geoplot(modis_polyshape);
        gp.EdgeAlpha = 0;
        % The color corresponds to the linear index_sort
        gp.FaceColor = C(ii==index_sort,:);
        gp.FaceAlpha = 0.7;

        

        hold on


    
end

% Plot the Vocals-Rex data
geoscatter(vocalsRex.latitude, vocalsRex.longitude, 10, "red",'*')


cb = colorbar;
%clim([min(modis.cloud.effRadius17(1:n_rows, 1:n_cols)), max(modis.cloud.effRadius17(1:n_rows, 1:n_cols))])
clim([min(modis.cloud.effRadius17(vocalsRex.modisIndex_minDist), [], 'all'),...
    max(modis.cloud.effRadius17(vocalsRex.modisIndex_minDist), [], 'all')])

set(get(cb, 'label'), 'string', '$r_e \; (\mu m)$','Interpreter','latex', 'Fontsize',28)
set(gca, 'FontSize',25)
set(gca, 'FontWeight', 'bold')
set(gcf, 'Position', [0 0 800 800])

% load the across and along pixel growth curve fits
if strcmp(whatComputer, 'anbu8374')==true

    title(['MODIS Effective Radius - ', modisFolder(97:end-1)],'Interpreter','latex', 'FontSize', 35)

elseif strcmp(whatComputer, 'andrewbuggee')==true

    title(['MODIS Effective Radius - ', modisFolder(113:end-1)],'Interpreter','latex', 'FontSize', 35)

end




end