function [f] = plot_MODIS_RGB(modis, inputs)


R = modis.EV.m250.radiance(:,:,2);      % 850 nm band is going to be our red hue
G = modis.EV.m250.radiance(:,:,1);      % 650 nm band is going to be our green hue
B = modis.EV.m500.radiance(:,:,3);      % 465 nm band is going to be our blue hue

numRows = size(R,1);
numCols = size(R,2);

% If there are any radiance values less than 0, set them to be 0
% R(R<0) = 0;
% G(G<0) = 0;
% B(B<0) = 0;

% Lets normalize each band and convert it to 8-bit image

R = uint8(R);
G = uint8(G);
B = uint8(B);

% create the iamge array
I = cat(3,R,G,B);

% for some reason we our aray is upside down. We want the equatorward
% direction to be south, and the polewrad direction to be north

I = rot90(rot90(I));

% lets make a plot!

f = figure; 
image(I); 
title(['True Color - ',inputs.INP_folderName(6:end-1)])
set(f, 'Position', [0 0 2.5*floor(numCols/10) 2*floor(numRows/10)])

end

