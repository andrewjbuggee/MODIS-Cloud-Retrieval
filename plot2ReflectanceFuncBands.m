%% ----- Plot two bands of the Reflectance Function against one another -----



% By Andrew J. Buggee

%%

function [] = plot2ReflectanceFuncBands(modis,R,inputs)

% extract inputs

re = inputs.re;
tau_c = inputs.tau_c;
pixel_row = inputs.pixel_row;
pixel_col = inputs.pixel_col;
bands2plot = inputs.bands2plot;

if length(bands2plot)~=2
    
    error('Can only plot two bands at a time')
    
end


figure;

% Plot values with constant particle radius, and vary the optical thickness
% the column dimension varies optical thickness

legendStr1 = cell(1,size(R,1));

for ii = 1:size(R,1)

    plot([R{ii,:,bands2plot(1)}],[R{ii,:,bands2plot(2)}]);
    
     hold on
          
     legendStr1{ii} = [num2str(re(ii)),' \mum'];


end

legendStr2 = cell(1,size(R,2));

text_x = zeros(1,size(R,2));
text_y = zeros(1,size(R,2));

% now plot the lines on constant optical thickness

for jj = 1:size(R,2)
    
    x = [R{:,jj,bands2plot(1)}];
    y = [R{:,jj,bands2plot(2)}];

    plot(x,y,'w--');
    
     hold on
     
     if jj<size(R,2)
        legendStr2{jj} = num2str(tau_c(jj));
     elseif jj == size(R,2)
         legendStr2{jj} = [num2str(tau_c(jj)),'   \tau_{c}'];
     end
     
     
    [~,idx] = min(y);
    text_x(jj) = x(idx)+ 0.005;
    text_y(jj) = y(idx) - 0.005;



end

    t = text(text_x, text_y, legendStr2, 'Color', 'w','fontsize',14,'fontweight','bold');

if bands2plot(1)<=2 && bands2plot(2)<=2
    
    xlabel(['Reflectance Function (',num2str(modis.EV.m250.bands.center(bands2plot(1))),' nm)']); 
    ylabel(['Reflectance Function (',num2str(modis.EV.m250.bands.center(bands2plot(2))),' nm)']);
    
elseif bands2plot(1)<=2 && bands2plot(2)>2
    
    yband = bands2plot(2) - 2;
    
    xlabel(['Reflectance Function (',num2str(modis.EV.m250.bands.center(bands2plot(1))),' nm)']); 
    ylabel(['Reflectance Function (',num2str(modis.EV.m500.bands.center(yband)),' nm)']);
    
    elseif bands2plot(1)>2 && bands2plot(2)>2
        
        xband = bands2plot(1) - 2;
        yband = bands2plot(2) - 2;
    
    xlabel(['Reflectance Function (',num2str(modis.EV.m500.bands.center(xband)),' nm)']); 
    ylabel(['Reflectance Function (',num2str(modis.EV.m500.bands.center(yband)),' nm)']);
    
    elseif bands2plot(1)>2 && bands2plot(2)<=2
        
        xband = bands2plot(1) - 2;
    
    xlabel(['Reflectance Function (',num2str(modis.EV.m500.bands.center(xband)),' nm)']); 
    ylabel(['Reflectance Function (',num2str(modis.EV.m250.bands.center(bands2plot(2))),' nm)']);
    
    
end


    grid on; grid minor
    legend(legendStr1,'location','best')
    
    title(['\theta_{0} = ',num2str(modis.solar.zenith(pixel_row,pixel_col)),'  \phi_{0} = ',...
        num2str(modis.solar.azimuth(pixel_row,pixel_col)),'  \theta = ',num2str(modis.sensor.zenith(pixel_row,pixel_col)),...
        '  \phi = ',num2str(modis.sensor.azimuth(pixel_row,pixel_col))]);
    
    
end


