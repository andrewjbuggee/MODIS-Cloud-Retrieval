%% ----- Plot two bands of the Reflectance Function against one another -----



% By Andrew J. Buggee

%%

function [] = plot2ReflectanceFuncBands(modis,R,modisInputs, pixels2use)

% extract inputs

re = modisInputs.re;
tau_c = modisInputs.tau_c;
pixel_row = pixels2use.res1km.row;
pixel_col = pixels2use.res1km.col;
bands2run = modisInputs.bands2run;
bands2plot = modisInputs.bands2plot;



if length(bands2plot)~=2

    error('Can only plot two bands at a time')

end


modis_band_vals = modisBands([modisInputs.bands2plot]);

% find the indices needed to plot the bands listed above
for ii = 1:length(bands2plot)
    index_bands(ii) = find(bands2run==bands2plot(ii));
end


% Create a legend string
legend_str = cell(1, length(tau_c) + 1 + length(re));
% set the first length(tau_c)+1 entries to be empty strings
for ss = 1:length(tau_c)+1
    legend_str{ss} = '';
end


% R is a 4-D matrix.
%   dim(R,1) = changing pixels
%   dim(R,2) = changing effective radius
%   dim(R,3) = changing optical thickness
%   dim(R,4) = changing spectral band



% Plot values with constant particle radius, and vary the optical thickness

% Plot across pixels
for pp = 1:size(R,1)

    figure;

    % Step through each effective radius. Each line represents the reflectance
    % at a constant droplet size, with varrying optical thickness.
    for rr = 1:length(re)

        plot(reshape(R(pp,rr,:,index_bands(1)),1,[]), reshape(R(pp,rr,:,index_bands(2)),1,[]),...
            '.-', 'MarkerSize',50,'LineWidth',1.5);
        hold on

        % store legend string for later
        legend_str{length(tau_c) + rr} = ['$r_e = $', num2str(re(rr))];

    end

    % set up color order for each curve
    colororder(mySavedColors(1:length(re),'fixed'));

    % Now plot the MODIS measurement on top
    plot(modis.EV1km.reflectance(pixel_row(pp), pixel_col(pp),modisInputs.bands2run(index_bands(1))),...
        modis.EV1km.reflectance(pixel_row(pp), pixel_col(pp),modisInputs.bands2run(index_bands(2))), 'diamond',...
        'MarkerSize',15, 'Color',mySavedColors(length(re)+1, 'fixed'), 'MarkerFaceColor','none');


    % ------ Plot lines of constant optical depth ------
    % Now step through each optical depth. Each line represents the reflectance
    % at a constant optical depth, with varrying effective radius.
    for tt = 1:length(tau_c)

        x = reshape(R(pp,:,tt,index_bands(1)),1,[]);
        y = reshape(R(pp,:,tt,index_bands(2)),1,[]);

        t = plot(x, y, 'LineStyle','--', 'Color','k');

        % add line label on plot
        if tt==1
            text(0.995*x(end), 0.8*y(end), num2str(tau_c(tt)),'Interpreter','latex',"FontSize",20, "FontWeight","bold")
            hold on
        else
            text(0.995*x(end), 0.8*y(end), num2str(tau_c(tt)),'Interpreter','latex',"FontSize",20, "FontWeight","bold")
            hold on
        end

        % place the dotted line below the lines of constant radius
        uistack(t, 'bottom');

    end

    % Add text to indicate the black lines are lines of constant optical
    % thickness
    text(1.05*x(end), 0.85*y(end), '$\tau_c$','Interpreter','latex',"FontSize",20, "FontWeight","bold")



    % set up plot stuff
    grid on; grid minor
    xlabel(['Reflectance ', num2str(modis_band_vals(1,1)), ' $nm$'],Interpreter='latex')
    ylabel(['Reflectance ', num2str(modis_band_vals(2,1)), ' $nm$'],Interpreter='latex')

    % set the last string entry to be MODIS value
    legend_str{end} = 'MODIS';

    legend(legend_str, 'Interpreter','latex','Location','best' , 'FontSize', 20, 'FontWeight','bold')
    title('Simulated Reflectance','Interpreter','latex')
    set(gcf,"Position", [0 0 1200 800])


end



end


