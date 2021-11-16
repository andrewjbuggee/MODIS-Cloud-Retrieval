%% ---- Get MODIS .INP Names based on pixel -----



% By Andrew J. Buggee

%%

function [inpNames] = getMODIS_INPnames_withClouds(solar,inputs)

% extract inputs

re = inputs.re;
tau_c = inputs.tau_c;
bands2run = inputs.bands2run;
pixel_row = inputs.pixel_row;
pixel_col = inputs.pixel_col;


rowName = num2str(pixel_row);
colName = num2str(pixel_col);

szaName = num2str(solar.zenith(pixel_row,pixel_col));
sazName = num2str(solar.azimuth(pixel_row,pixel_col));

inpNames = cell(length(re),length(tau_c),length(bands2run));

for ii = 1:length(bands2run)
    
    for jj = 1:length(re)
        
        for kk = 1:length(tau_c)
    
        bandName = num2str(bands2run(ii));
        inpNames{jj,kk,ii} = ['pixel_',rowName,'r_',colName,...
                 'c_sza_',szaName,'_saz_',sazName,'_band_',bandName,...
                 '_r_',num2str(re(jj)),'_T_',num2str(tau_c(kk)),'.INP'];
    
        end

    end
    
end





end
