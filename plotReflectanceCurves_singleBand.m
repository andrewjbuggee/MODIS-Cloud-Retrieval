%% Plot reflectance curves of constant droplet radius


% By Andrew J. Buggee
%%

function f = plotReflectanceCurves_singleBand(R,inputs,pixels2use)

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
        
        f(pp) = figure;
        
        for ii = 1:size(bands2plot,1)
            
            
            
            for jj = 1:size(bands2plot,2)
                
                % find indices for bands 2 plot
                bands2plot_index = bands2plot(ii,jj) == bands2run;
                
                subplot(1,size(bands2plot,2),jj)
                
                bandVals = modisBands(bands2plot(ii,jj));
                
                lgnd_str = cell(1,length(re));
                
                for kk = 1:length(re)
                    
                    
                    reflectance = R(pp,kk,:,bands2plot_index);
                    plot(tau_c,reflectance(:));
                    hold on
                    lgnd_str{kk} = ['r_{e} = ',num2str(re(kk)),' \mum'];
                    set(gcf, 'Position', [0 0 1500 600])

                end
                
                
                title([num2str(bandVals(ii,jj)),' nm'])
                xlabel('\tau_{c}'); ylabel('Reflectance')
                grid on; grid minor
                
                if jj == 1
                    legend(lgnd_str,'Location','best')
                    
                    dim = [.5 0 .3 .3];
                    str = ['sza = ',num2str(sza(pp)),' saz = ',num2str(saz(pp)),' vza = ',num2str(vza(pp)),...
                        ' vaz = ',num2str(vaz(pp))];
                    annotation('textbox',dim,'String',str,'FitBoxToText','on','Color','k',...
                        'FontWeight','bold','FontSize',14, 'EdgeColor', 'w');
                end
                
                
            end
            
        end
        
        
    end
    
    
elseif num_pixels > 3
    
    % if there are a bunch of pixels, we will just grab a random subset of
    % 3 to plot
    rand_index = randsample(num_pixels,3); % random sampling without replacement
    
    for pp = 1:length(rand_index)
        
        for ii = 1:size(bands2plot,1)
            
            f(pp) = figure;
            
            for jj = 1:size(bands2plot,2)
                
                 % find indices for bands 2 plot
                bands2plot_index = bands2plot(ii,jj) == bands2run;
                
                subplot(1,size(bands2plot,2),jj)
                
                bandVals = modisBands(bands2plot(ii,jj));
                
                lgnd_str = cell(1,length(re));
                
                for kk = 1:length(re)
                    
                    reflectance = R(rand_index(pp),kk,:,bands2plot_index);
                    plot(tau_c,reflectance(:));
                    hold on
                    lgnd_str{kk} = ['r_{e} = ',num2str(re(kk)),' \mum'];
                    set(gcf, 'Position', [0 0 1500 600])

                    
                end
                
                
                title([num2str(bandVals(1)),' nm'])
                xlabel('\tau_{c}'); ylabel('Reflectance')
                grid on; grid minor
                
                if jj == 1
                    legend(lgnd_str,'Location','best')
                    
                    dim = [.5 0 .3 .3];
                    str = ['sza = ',num2str(sza(rand_index(pp))),' saz = ',num2str(saz(rand_index(pp))),' vza = ',num2str(vza(rand_index(pp))),...
                        ' vaz = ',num2str(vaz(rand_index(pp)))];
                    annotation('textbox',dim,'String',str,'FitBoxToText','on','Color','k',...
                        'FontWeight','bold','FontSize',14, 'EdgeColor','w');
                end
                
                
            end
            
        end
        
    end
    
    
    
    
    
end




