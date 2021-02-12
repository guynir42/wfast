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
%   -weights: give the relative weight foreach measurement, or a single
%             value for all measurements. Will be used as weights for the
%             fit and also for scaling the chi2 result. 
%             The weights are treated as the inverse errors, and are
%             ignored if the errors are given.
%   -errors: give the standard deviation for each measurement. 
%            This is used to scale the chi2 and weigh the fit. 
%            This input overrides the "weights" input. 
%   -order: highest polynomial order in both x and y. Default is 2. 
%   -iterations: how many iterations (to remove outliers). Default is 1. 
%   -sigma: remove outliers above this threshold. Default is 5. 
%   -x_values and y_values: give range or list of values for making the
%                           surface map output and for plotting. 
%                           Default is min/max of x and y values. 
%   -size: the number of pixels in the surface map output. Default 1000. 
%   -double: force the data to be in double-precision. Default is true. 
%   -plot: show the fit and the data points. Default is false. 
%   -pause: when plotting, how many seconds to wait between iterations.
%           Default is one second. 
%   -axes: the graphic axes to plot into. Default is gca(). 
%   -font_size: the axes font size when plotting. 
%   -marker_size: the size of the little circles for the measurements, on 
%                 the scatter plot (when plot=1). Default is 40. 
%   -scale: set the contrast limits (CLim) or scale of the plot. Default is
%           to set the limits automatically. 
%   -warnings: When true, will show Rand Deficient Matrix warnings. 
%              Default is false (warnings are ignored). 



    if nargin<3, help('util.fit.surf_poly'); return; end

    input = util.text.InputVars;
    input.use_ordered_numeric = 1; 
    input.input_var('weights', []); 
    input.input_var('errors', []); 
    input.input_var('order', 2); 
    input.input_var('iterations', 1); 
    input.input_var('sigma', 5); 
    input.input_var('x_values', []); 
    input.input_var('y_values', []); 
    input.input_var('size', 1000, 'im_size', 'image_size'); 
    input.input_var('double', true); 
    input.input_var('plot', false); 
    input.input_var('pause', 1, 'duration'); 
    input.input_var('ax', [], 'axes', 'axis'); 
    input.input_var('font_size', 18); 
    input.input_var('marker_size', 40); 
    input.input_var('scale', [], 'CLim'); 
    input.input_var('warnings', false); % if turned on, will show warnings like Rand Deficient Matrix
    input.scan_vars(varargin{:}); 
    
    
    warn = warning('query', 'MATLAB:lscov:RankDefDesignMat'); % check the current state of the warning
    prev_state = warn.state; 
    new_state = 'off'; 
    if input.warnings, new_state = 'on'; end
    warning(new_state, 'MATLAB:lscov:RankDefDesignMat')
    on_cleanup = onCleanup(@() warning(prev_state, 'MATLAB:lscov:RankDefDesignMat')); 
    
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
    
    w = [];
    
    if ~isempty(input.weights) 
        w = input.weights;
    end
    
    if ~isempty(input.errors)
        w = 1./input.errors;
    end
    
    if ~isempty(w)
        
        if isscalar(w)
            w = ones(size(x)).*w;
        end
        
        % maybe later we will add a case where the same weights are replicated for all columns in x, y, and v
        if ~isequal(size(x), size(w))
            error('Must input weights with the same size as x,y, and v!');
        end
        
        if isvector(w)
            w = util.vec.tocolumn(w); 
        end
        
        w(w<0) = 0; 
        
    end
    
    if input.double
        x = double(x); 
        y = double(y); 
        v = double(v); 
        w = double(w); 
    end
    
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
    
    N = length(coeff_names);
    
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
        
        X(bad_indices) = [];
        Y(bad_indices) = [];
        V(bad_indices) = [];

        if ~isempty(w)
            W = w(:,ii);                    
            W(bad_indices) = [];
        end
        
        for jj = 1:input.iterations

            if length(nan_indices) + N >= size(x,1)
                new_result = make_empty_result(coeff_names, coeff_powers_x, coeff_powers_y, x(:,ii), y(:,ii), v(:,ii), w(:,ii));
                break;
            end
            
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
            new_result.coeff_powers_x = coeff_powers_x;
            new_result.coeff_powers_y = coeff_powers_y;
            new_result.coeffs = coeffs;

            new_result.model = ''; 
            for kk = 1:length(coeff_names)
                new_result.model = sprintf('%s%+f', new_result.model, coeffs(kk)); 
                if ~isempty(coeff_names{kk})
                    new_result.model = sprintf('%s.*%s', new_result.model, coeff_names{kk});
                end
            end

            new_result.func = str2func(['@(x,y) ' new_result.model]); 

            new_result.x = x(:,ii);
            new_result.y = y(:,ii);
            new_result.v = v(:,ii);
            if ~isempty(w)
                new_result.w = w(:,ii);
            end
            
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
                new_result.chi2 = nansum((reduced_residuals.*W).^2); 
            end
            
            new_result.ndof = length(V) - N; 
            
            new_result.var = new_result.chi2./new_result.ndof; % the measurement error estimate
            
            if input.plot

                I = make_image(coeffs, coeff_powers_x, coeff_powers_y, input.x_values, input.y_values);
                
                imagesc(input.ax, input.x_values, input.y_values, I); 
                
                hold(input.ax, 'on'); 
                
                if isempty(W)
                    marker_sizes = input.marker_size;
                else
                    marker_sizes = input.marker_size.*w(:,ii)./nanmean(w(:,ii)) + 1;
                end
                
                scatter(input.ax, x(:,ii), y(:,ii), marker_sizes, v(:,ii), 'o', 'filled','MarkerEdgeColor','k'); 
                scatter(input.ax, x(bad_indices,ii), y(bad_indices,ii), input.marker_size, 'r', 'x'); 
                hold(input.ax, 'off'); 
                
                colorbar(input.ax); 
                input.ax.FontSize = input.font_size;
                
                if ~isempty(input.scale)
                    input.ax.CLim = input.scale; 
                end
                
                pause(input.pause); 
            
            end
            
            if isempty(w)
                idx = find(abs(new_result.residuals) > input.sigma.*sqrt(new_result.var));
            else
                idx = find(abs(new_result.residuals.*w(:,ii)) > input.sigma.*sqrt(new_result.var));
            end
            
            idx = unique([nan_indices; idx]); % add the measurements where we found NaN values
            
            if isequal(idx, bad_indices) || length(idx)==size(x,1)
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

function new_result = make_empty_result(coeff_names, coeff_powers_x, coeff_powers_y,x,y,v,w)
    
    coeffs = NaN(length(coeff_names),1); 
    new_result.coeff_names = coeff_names;
    new_result.coeff_powers_x = coeff_powers_x;
    new_result.coeff_powers_y = coeff_powers_y;
    new_result.coeffs = coeffs;
    new_result.model = ''; 
    for kk = 1:length(coeff_names)
        new_result.model = sprintf('%s%+f', new_result.model, coeffs(kk)); 
        if ~isempty(coeff_names{kk})
            new_result.model = sprintf('%s.*%s', new_result.model, coeff_names{kk});
        end
    end

    new_result.func = @(x,y) NaN; 

    new_result.x = x;
    new_result.y = y;
    new_result.v = v;
    if ~isempty(w)
        new_result.w = w;
    end
    new_result.vm = NaN(size(v,1),1); 

    new_result.residuals = NaN(size(v,1),1);
    new_result.chi2 = NaN;
    new_result.ndof = length(v) - length(coeffs); 
    new_result.var = NaN;
    
end





