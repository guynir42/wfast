function [results, surf_image] = surf_poly(x,y,v,varargin)
% Usage: [results, surf_image] = surf_poly(x,y,v,varargin)
% Fit a polynomial surface in x and y to the values in v. 
% 
% INPUTS: give the x, y and value at each measured point. 
%
% OUTPUTS: -results: a struct containing the ocefficients etc. 
%          -surf_image: a map of the surface using the fit results. 
%
% OPTIONAL PARAMETERS:
%   -weights: give the error estimate for each measurement, or a single
%             value for all measurements. Will be used as weights for the
%             fit and also for scaling the chi2 result. 
%   -order: highest polynomial order in both x and y. Default is 2. 
%   -iterations: how many iterations (to remove outliers). Default is 1. 
%   -sigma: remove outliers above this threshold. Default is 5. 
%   -x_values and y_values: give range or list of values for making the
%                           surface map output and for plotting. 
%                           Default is min/max of x and y values. 
%   -size: the number of pixels in the surface map output. Default 1000. 
%   -plot: show the fit and the data points. Default is false. 
%   -pause: when plotting, how many seconds to wait between iterations.
%           Default is one second. 
%   -axes: the graphic axes to plot into. Default is gca(). 
%   -font_size: the axes font size when plotting. 
%           



    if nargin<3, help('util.fit.surf_poly'); return; end

    input = util.text.InputVars;
    input.use_ordered_numeric = 1; 
    input.input_var('weights', [], 'errors'); 
    input.input_var('order', 2); 
    input.input_var('iterations', 1); 
    input.input_var('sigma', 5); 
    input.input_var('x_values', []); 
    input.input_var('y_values', []); 
    input.input_var('size', 1000, 'im_size', 'image_size'); 
    input.input_var('plot', false); 
    input.input_var('pause', 1, 'duration'); 
    input.input_var('ax', [], 'axes', 'axis'); 
    input.input_var('font_size', 18); 
    input.scan_vars(varargin{:}); 
    
    if input.plot && isempty(input.ax)
        input.ax = gca;
    end
    
    if input.plot || nargout>1 % need to produce an image map
        
        s = util.vec.imsize(input.size); % make sure it is a 2 element vector
        
        if isempty(input.x_values)
            input.x_values = linspace(nanmin(x(:)), nanmax(x(:)), s(2)); 
        elseif length(input.x_values)==2
            input.x_values = linspace(input.x_values(1), input.x_values(2), s(2)); 
        end
        
        if isempty(input.y_values)
            input.y_values = linspace(nanmin(y(:)), nanmax(y(:)), s(1)); 
        elseif length(input.y_values)==2
            input.y_values = linspace(input.y_values(1), input.y_values(2), s(1)); 
        end
        
    end
    
    if ~isequal(size(x), size(y)) || ~isequal(size(x), size(v))
        error('Must input x,y, and v of the same size!'); 
    end
    
    if isvector(x)
        x = util.vec.tocolumn(x); 
        y = util.vec.tocolumn(y); 
        v = util.vec.tocolumn(v);
    end
    
    if ~isempty(input.weights) 
        
        if isscalar(input.weights)
            input.weights = ones(size(x)).*input.weights;
        end
        
        % maybe later we will add a case where the same weights are replicated for all columns in x, y, and v
        if ~isequal(size(x), size(input.weights))
            error('Must input weights with the same size as x,y, and v!');
        end
        
        if isvector(input.weights)
            input.weights = util.vec.tocolumn(input.weights); 
        end
        
    end
    
    w = input.weights; 
    
    S = size(x); % keep track of the original size (unless given a row vector, but then it doesn't matter)
    
    % linearize all higher dimensions
    x = x(:,:); 
    y = y(:,:); 
    v = v(:,:); 
    w = w(:,:); 
    
    coeff_names = {''};
    coeff_powers_x = 0;
    coeff_powers_y = 0;
    
    for ord = 1:input.order
        
        for ii = 0:ord
            
            coeff_powers_x(end+1) = ord - ii;
            coeff_powers_y(end+1) = ii; 
            
            if coeff_powers_x(end)==0
                name_x = '';
            elseif coeff_powers_x(end)==1
                name_x = 'x';
            else
                name_x = sprintf('x.^%d', coeff_powers_x(end));
            end
            
            if coeff_powers_y(end)==0
                name_y = '';
            elseif coeff_powers_y(end)==1
                name_y = 'y';
            else
                name_y = sprintf('y.^%d', coeff_powers_y(end));
            end
            
            if ~isempty(name_x) && ~isempty(name_y)
                name = [name_x '.*' name_y]; 
            elseif ~isempty(name_x)
                name = name_x;
            elseif ~isempty(name_y)
                name = name_y;
            else
                name = '';
            end
            
            coeff_names{end+1} = name;
            
        end
        
    end
    
    for ii = 1:size(x,2)
        
        % we can exclude some of these
        X = x(:,ii);
        Y = y(:,ii);
        V = v(:,ii);        
        if ~isempty(w)
            W = w(:,ii);
        else
            W = [];
        end
        
        nan_indices = isnan(x(:,ii)) | isnan(y(:,ii)) | isnan(v(:,ii)); 
        
        if ~isempty(W)
            nan_indices = nan_indices | isnan(w(:,ii)); 
        end
        
        nan_indices = find(nan_indices); 
        
        bad_indices = nan_indices; % indices of outlier measurements or NaNs
        
        for jj = 1:input.iterations

            A = []; % the design matrix is built up from powers of x and y
            for kk = 1:length(coeff_powers_x)
                A = [A, X.^coeff_powers_x(kk).*Y.^coeff_powers_y(kk)]; 
            end

            if isempty(input.weights)
                coeffs = lscov(A,V); 
            else
                coeffs = lscov(A,V,W); 
            end

            new_result.coeff_names = coeff_names;
            new_result.coeffs = coeffs;

            new_result.model = ''; 
            for kk = 1:length(coeff_names)
                new_result.model = sprintf('%s%+f', new_result.model, coeffs(kk)); 
                if ~isempty(coeff_names{kk})
                    new_result.model = sprintf('%s.*%s', new_result.model, coeff_names{kk});
                end
            end

            new_result.func = str2func(['@(x,y) ' new_result.model]); 

            new_result.x = x;
            new_result.y = y;
            new_result.v = v;
            new_result.w = w;
            
            % produce the values at all x,y points based on the model
            B = 0;
            for kk = 1:length(coeff_powers_x)
                B = B + coeffs(kk).*x(:,ii).^coeff_powers_x(kk).*y(:,ii).^coeff_powers_y(kk); % these results include the outlier positions
            end       
            
            new_result.vm = B;
            new_result.residuals = v(:,ii)-B; 
            
            B2 = B;
            B2(bad_indices) = []; % remove points where outliers were found
            
            if isempty(bad_indices)
                reduced_residuals = new_result.residuals; % just use all points
            else
                reduced_residuals = V - B2; 
            end
            
            if isempty(W)
                new_result.chi2 = nansum(reduced_residuals.^2); 
            else
                new_result.chi2 = nansum((reduced_residuals./W).^2); 
            end
            
            new_result.ndof = length(V) - length(coeffs); 
            
            new_result.var = new_result.chi2./new_result.ndof; % the measurement error estimate
            
            if input.plot

                I = make_image(coeffs, coeff_powers_x, coeff_powers_y, input.x_values, input.y_values);
                
                imagesc(input.ax, input.x_values, input.y_values, I); 
                
                hold(input.ax, 'on'); 
                
                if isempty(W)
                    marker_sizes = 25;
                else
                    marker_sizes = 10.*w(:,ii)./nanmean(w(:,ii)); 
                end
                
                scatter(input.ax, x(:,ii), y(:,ii), marker_sizes, v(:,ii), 'o', 'filled'); 
                scatter(input.ax, x(bad_indices,ii), y(bad_indices,ii), 25, 'r', 'x'); 
                hold(input.ax, 'off'); 
                
                colorbar(input.ax); 
                input.ax.FontSize = input.font_size;
                
                pause(input.pause); 
            
            end
            
            idx = find(abs(new_result.residuals./w(:,ii)) > input.sigma.*sqrt(new_result.var));
            idx = unique([nan_indices; idx]); % add the measurements where we found NaN values
            
            if isequal(idx, bad_indices)
                break;
            else
                
                X = x(:,ii);
                Y = y(:,ii);
                V = v(:,ii);        
                
                X(idx) = [];
                Y(idx) = [];
                V(idx) = [];
                
                if ~isempty(w)
                    W = w(:,ii);                    
                    W(idx) = [];
                end
                
                bad_indices = idx; 
                
            end
            
            new_result.outlier_indices = bad_indices;

        end % for jj (iterations)
        
        results(ii) = new_result;
        
    end % for ii (columns/value sets)
    
    if nargout>1 % create an image from the surface fit
        
        surf_image = zeros(length(input.x_values), length(input.x_values), length(results)); 
        
        for ii = 1:size(x,2)
            surf_image(:,:,ii) = make_image(results(ii).coeffs, coeff_powers_x, coeff_powers_y, input.x_values, input.y_values); 
        end
        
        surf_image = reshape(surf_image, [size(surf_image,1), size(surf_image,2), S(2:end)]); 
        
    end
    
    results = reshape(results, [1, S(2:end)]); 

end

function I = make_image(coeff_values, coeffs_x, coeffs_y, x_values, y_values)
    
    [X,Y] = meshgrid(x_values, y_values); 

    I = 0;
    for kk = 1:length(coeffs_x)
        I = I + coeff_values(kk).*X.^coeffs_x(kk).*Y.^coeffs_y(kk);
    end         

end







