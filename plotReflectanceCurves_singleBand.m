%% Plot reflectance curves of constant droplet radius


% By Andrew J. Buggee
%%

function [] = plotReflectanceCurves_singleBand(R,inputs,pixels2use)

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



if num_pixels <=3
    
    for pp = 1:num_pixels
        
        figure;
        
        for ii = 1:size(bands2plot,1)
            
            
            
            for jj = 1:size(bands2plot,2)
                
                % find indices for bands 2 plot
                bands2plot_index = bands2plot(ii,jj) == bands2run;
                
                subplot(1,size(bands2plot,2),jj)
                
                band_center = modisBandsCenter(bands2plot(ii,jj));
                
                lgnd_str = cell(1,length(re));
                
                for kk = 1:length(re)
                    
                    
                    reflectance = R(pp,kk,:,bands2plot_index(ii,:));
                    plot(tau_c,reflectance(:));
                    hold on
                    lgnd_str{kk} = ['r_{e} = ',num2str(re(kk)),' \mum'];
                    
                end
                
                
                title(['Reflectance - ',num2str(band_center),' \mum'])
                xlabel('\tau_{c}'); ylabel('r_{e} (\mum)')
                grid on; grid minor
                
                if jj == 1
                    legend(lgnd_str,'Location','best')
                    
                    dim = [.5 0 .3 .3];
                    str = ['sza = ',num2str(sza(pp)),' saz = ',num2str(saz(pp)),' vza = ',num2str(vza(pp)),...
                        ' vaz = ',num2str(vaz(pp))];
                    annotation('textbox',dim,'String',str,'FitBoxToText','on','Color','white',...
                        'FontWeight','bold','FontSize',14);
                end
                
                
            end
            
        end
        
        
    end
    
    
elseif num_pixels > 3
    
    % if there are a bunch of pixels, we will just grab a random subset of
    % 3 to plot
    rand_index = randi(num_pixels,1,3);
    
    for pp = 1:num_pixels
        
        for ii = 1:size(bands2plot,1)
            
            figure;
            
            for jj = 1:size(bands2plot,2)
                
                subplot(1,size(bands2plot,2),jj)
                
                band_center = modisBandsCenter(bands2plot(ii,jj));
                
                lgnd_str = cell(1,length(re));
                
                for kk = 1:length(re)
                    
                    plot(tau_c,[R(rand_index(pp),kk,:,bands2plot(jj))]);
                    hold on
                    lgnd_str{kk} = ['r_{e} = ',num2str(re(kk)),' \mum'];
                    
                end
                
                
                title(['Reflectance - ',num2str(band_center),' \mum'])
                xlabel('\tau_{c}'); ylabel('r_{e} (\mum)')
                grid on; grid minor
                
                if jj == 1
                    legend(lgnd_str,'Location','best')
                    
                    dim = [.5 0 .3 .3];
                    str = ['sza = ',num2str(sza(pp)),' saz = ',num2str(saz(pp)),' vza = ',num2str(vza(pp)),...
                        ' vaz = ',num2str(vaz(pp))];
                    annotation('textbox',dim,'String',str,'FitBoxToText','on','Color','white',...
                        'FontWeight','bold','FontSize',14);
                end
                
                
            end
            
        end
        
    end
    
    
    
    
    
end




