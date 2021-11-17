%% Create Reflectance Contour Plots



% By Andrew J. Buggee
%%

function [] = plotReflectanceContours(R,inputs,pixels2use)

bands2run = inputs.bands2run;
bands2plot = inputs.bands2plot;
re = inputs.re;
tau_c = inputs.tau_c;
num_pixels = inputs.pixels.num_2calculate;

% extract pixel geometry
sza = pixels2use.res1km.geometry.sza; % solar zenith angle
saz = pixels2use.res1km.geometry.saz; % solar azimuth angle
vza = acosd(pixels2use.res1km.geometry.umu); % viewing zenith angle
vaz = pixels2use.res1km.geometry.phi; % viewing azimuth angle

for pp = 1:num_pixels
    
    figure;
    
    for ii = 1:size(bands2plot,1)
        
        
        
        for jj = 1:size(bands2plot,2)
            
            % find indices for bands 2 plot
            bands2plot_index = bands2plot(ii,jj) == bands2run;
            
            band_center = modisBandsCenter(bands2plot(ii,jj));
            subplot(1,size(bands2plot,2),jj)
            
            reflectance = reshape(R(pp,:,:,bands2plot_index),length(re),length(tau_c));
            contourf(tau_c,re,reflectance,'ShowText','on'); colorbar;
            title(['Reflectance - ',num2str(band_center),' \mum'])
            xlabel('\tau_{c}'); ylabel('r_{e} (\mum)')
            
        end
        
        
    end
    
    dim = [.5 0 .3 .3];
    str = ['sza = ',num2str(sza(pp)),' saz = ',num2str(saz(pp)),' vza = ',num2str(vza(pp)),...
        ' vaz = ',num2str(vaz(pp))];
    annotation('textbox',dim,'String',str,'FitBoxToText','on','Color','white',...
        'FontWeight','bold','FontSize',14);
    
    
end









end