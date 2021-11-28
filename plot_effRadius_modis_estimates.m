%% ----- Compare the effectuve radius calculate by modis, and my algorithm -----



% By Andrew J. Buggee
%%

function [] = plot_effRadius_modis_estimates(truth_estimate_table)

% extract the modis estimate and my calculation estimates
modis_R17 = truth_estimate_table.modisR17;
est_R17 = truth_estimate_table.estR17;

abs_diff = truth_estimate_table.absDiffR; % the absolute difference between my estimate and the modis estimate
avg_abs_diff = mean(abs_diff);



% find the minimum and maximum values to create a y=x line

min_est = min(est_R17);
min_modis = min(modis_R17);

max_est = max(est_R17);
max_modis = max(modis_R17);

min_global = min([min_est,min_modis]);

max_global = min([max_est,max_modis]);

x = linspace((0.9 * min_global),(1.1*max_global),150);


figure; plot(x,x,'w-','Linewidth',1)
hold on; grid on; grid minor
plot(est_R17,modis_R17,'m.')
xlabel('My estimate - r_{e} (\mum)')
ylabel('MODIS estimate - r_{e} (\mum)') 


% find the indices of estiamtes that are furthest from their modis
% counterpart

num2find = 10;
index_2find = zeros(1,num2find);

for ii = 1:num2find
    
    [~, index_2find(ii)] = max(abs_diff);
    
    abs_diff(index_2find(ii)) = 0; % set it to a value that will never be chosen!
    
end
    

% what do I want to look at with the estimates that deviate the most from
% the modis values? 

% first lets plot all pixels again, and then highlight the 10 that deviate
% the most. This allows us to see how they stack up with the full set

figure; plot(x,x,'w-','Linewidth',1)
hold on; grid on; grid minor
plot(est_R17,modis_R17,'m.')
xlabel('My estimate - r_{e} (\mum)')
ylabel('MODIS estimate - r_{e} (\mum)') 
plot(est_R17(index_2find),modis_R17(index_2find),'c.')
legend('Perfect Fit','all pixels',[num2str(num2find),' furthest from line'],'Location','best')
title(['Mean Abs Difference: ',num2str(avg_abs_diff),' \mum']) 



% find and remove values of tau that modis deems to be greater than 80.
% This is the upper limit I set in my look up tables. 

index_80 = truth_estimate_table.modisT17>80;
est_R17_80 = est_R17(~index_80);
modis_R17_80 = modis_R17(~index_80);

% redefine the absolute difference since we altered some values up above.
% Also, ignore any values where the tau is greater than 80

abs_diff = truth_estimate_table.absDiffR(~index_80); % the absolute difference between my estimate and the modis estimate
num2find = 10;
index_2find = zeros(1,num2find);

for ii = 1:num2find
    
    [~, index_2find(ii)] = max(abs_diff);
    
    abs_diff(index_2find(ii)) = 0; % set it to a value that will never be chosen!
    
end

figure; plot(x,x,'w-','Linewidth',1)
hold on; grid on; grid minor
plot(est_R17_80,modis_R17_80,'m.')
xlabel('My estimate - r_{e} (\mum)')
ylabel('MODIS estimate - r_{e} (\mum)')  
plot(est_R17_80(index_2find),modis_R17_80(index_2find),'c.')
legend('Perfect Fit','all pixels',[num2str(num2find),' furthest from line'],'Location','best')
title(['Mean Abs Difference: ',num2str(mean(abs_diff))]) 



end