%% ----- Reflectance Function - Coarse Grid Search -----


% right now, this can only use two wavelengths to perform a minimum least
% squares grid search

% By Andrew J. Buggee

%%

function [minVals] = leastSquaresGridSearch(modisRefl,modelData,inputs,pixels2use)



% extract inputs
re = inputs.re;
tau_c = inputs.tau_c;
interpGridScaleFactor = inputs.interpGridScaleFactor;
bands2search = inputs.bands2search; % bands to determine optical properties
bands2run = inputs.bands2run; % bands to run through uvspec

% save calculations
saveCalcs_filename = inputs.saveCalculations_fileName;

index = zeros(size(bands2search,1),length(bands2run));

for ii = 1:size(bands2search,1)
    for jj = 1:size(bands2search,2)
        index_hold = bands2search(ii,jj) == inputs.bands2run;
        index(ii,:) = index(ii,:) + index_hold;
    end
end

band_index = logical(index);

numPixels = size(modelData,1); % number of pixels to preform grid search on



% now extract the mid point of each MODIS band being modeled
bandVals = modisBands(bands2search);

% lets interpolate the model data to increase our grid

% in our model data, the column space, whixh spans the x direction, varies
% with tau. The row space, which spans the y direction, vaires with re

[T0,Re0] = meshgrid(tau_c,re);

% now lets expand the two variables in question to the preferred size

newTau_c = linspace(min(tau_c),max(tau_c),interpGridScaleFactor*length(tau_c));
new_re = linspace(min(re),max(re),interpGridScaleFactor*length(re));

% set up the new grid to interpolate on
[T,Re] = meshgrid(newTau_c,new_re);

% SWITCH TO 1KM DATA
% but for now, lets propose a simple algebraic fix so taht we select
% the proper pixel in our 500 meter data set
pixelIndex_500m = 1:4:((numPixels-1)*4 +1);


for pp = 1:numPixels
    
    
    
    
    % lets shape our data in an easy to use format
    % first extract the data from the bands of interest
    
    % lets compare with non-interpolated data and use the 1km resolution
    % modis reflected data, right now, is in 500 meter resolution. So we
    % have to use the 500 meter pixel indexes
    for bb = 1:size(bands2search,1)
        
        observations = modisRefl(pixels2use.res1km.row(pp),pixels2use.res1km.col(pp),band_index(bb,:));
        
        if iscell(modelData) == true
            
            modelData_vec = [modelData{pp,:,:,band_index(bb,:)}];
            modelData_vec = reshape(modelData_vec,size(modelData,2),size(modelData,3),length(bands2search));
            
            
            
        elseif isnumeric(modelData) == true
            
            modelData_vec = modelData(pp,:,:,band_index(bb,:));
            modelData_vec = reshape(modelData_vec,size(modelData,2),size(modelData,3),size(bands2search,2));
            
            
        end
        
        % preform 2D interpolation
        newModelData = zeros(length(new_re),length(newTau_c),size(modelData_vec,3));
        
        for ii = 1:size(modelData_vec,3)
            newModelData(:,:,ii) = interp2(T0,Re0,modelData_vec(:,:,ii),T,Re);
        end
        
        % finally, lets now rescale the observation to be on a flat plane with the
        % same dimensions as the interpolated model data
        
        observations_newGrid = repmat(observations,length(new_re),length(newTau_c));
        
        
        %% ---- lets view the surfaces of the model -----
        
        band2Plot = 2;
        
        if inputs.flags.plotMLS_figures == true
            surfPlots4modisModel_andObs(T,Re,newModelData(:,:,band2Plot),observations_newGrid(:,:,band2Plot),bandVals(band2Plot,1))
        end
        %% ----- Least Squares Difference ------
        
        % take the difference between the model and the observaions at each
        % wavelength. The we square the difference and sum along wavelenghts. Then
        % we take the sqrt
        
        leastSquaresGrid = sqrt(sum((newModelData - observations_newGrid).^2,3));
        
        % If we assume that there is a unique global minimum, we can search for it
        % using the minimum function. If we have more complicated scenario we need
        % to use the optimization tools to search for global and local minima
        
        [minVals.minLSD(bb,pp),index] = min(leastSquaresGrid,[],'all','linear');
        [row,col] = ind2sub(size(leastSquaresGrid),index);
        
        minVals.minR(bb,pp) = Re(row,col);
        minVals.minT(bb,pp) = T(row,col);
        
        % lets look at the least squares grid
        if inputs.flags.plotMLS_figures == true
            
            figure; s = surf(T,Re,leastSquaresGrid);
            xlabel('Optical Depth')
            ylabel('Effective Radius (\mum)')
            zlabel(['Least Squares Difference'])
            s.EdgeColor = 'k';
            s.EdgeAlpha = 0.2;
            colorbar
            title('Error Function')
        end
        
    end
    
    
    save(saveCalcs_filename,"minVals",'-append'); % save inputSettings to the same folder as the input and output file
    
    
    
    %% ---- Global Minimum Search -----
    
    % -- PROBLEM -- I need to create a proper function handle, rather than
    % using a surface fit object to create a function handle. This is giving me
    % erroneous resutls when using the global search function.
    
    % vec_T = reshape(T,[],1);
    % vec_Re = reshape(Re,[],1);
    % vec_leastSquares = reshape(leastSquaresGrid,[],1);
    %
    % leastSquares_fit = fit([vec_T,vec_Re],vec_leastSquares,'poly23');
    %
    %
    % % create a function handle using the surface fit
    %
    % leastSquares_fh = @(T,r) leastSquares_fit(T,r);
    %
    % % lets define a starting point for the minimization solver
    % x0 = [0,0]; % start at 0 r and 0 tau
    % % create lower and upper bounds
    % lb = [0,0];
    % ub = [40,80];
    %
    % % define the problem to solve
    %
    %
    % problem = createOptimProblem('fmincon','objective',leastSquares_fh,'x0',x0,'lb',lb,'ub',ub,...
    %                             'options',optimoptions(@fmincon,'Algorithm','sqp','Display','off'));
    %
    % % find local minimums
    % gs = GlobalSearch('Display','iter');
    %
    % [xmin,fmin] = run(gs,problem)
    
    
end






