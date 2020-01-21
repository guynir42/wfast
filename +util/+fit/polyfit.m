function results = polyfit(x,y,varargin)
% Usage: results = polyfit(x,y,varagin)
% Fits the data in x and y to a polynomial, using linear least squares and 
% removing outliers iteratively. 
% 
% Inputs: x and y, e.g., x is the time and y the flux. 
% OPTIONAL ARGUMENTS: 
%   -order: the highest power of the polynomial (default 2, e.g., ax^2+bx+c)
%   -sigma: how many times the noise rms is considered outlier (default 3). 
%   -iterations: how many times to refit after removing outliers (default 3). 
%   -variances: if you have the error^2 per sample, it can be given to the fitter
%   -double: use double (64 bit) precision instead of single (32 bit). Default is false. 
%   -plotting: make the fitter plot the results of each iteration for debugging. 
%   -axes: what axes object to plot to. 
%   -pause: length of pause after each plot. Default 0.3 seconds. 

    if nargin==0, help('util.fit.polyfit'); return; end

    input = util.text.InputVars;
    input.input_var('order', 2); % polynomial order
    input.input_var('sigma', 3); % above this threshold is considered outliers
    input.input_var('iterations', 3); % how many iterations of removing outliers
    input.input_var('variances', []); % assume no covariance, only different variance per sample
    input.input_var('double', false); % use double precision
    input.input_var('plotting', 0);
    input.input_var('axes', [], 'axis');
    input.input_var('pause', 0.3);
    % add other options
    input.scan_vars(varargin{:});
    
    if input.double
        y = double(y);
    end
    
    if ndims(x)>2
        error('Cannot handle more than 2D for x');
    end
    
    if ndims(y)>2
        error('Cannot handle more than 2D for y');
    end
    
    if isvector(x)
        x = util.vec.tocolumn(x);
    end
    
    if isvector(y)
        y = util.vec.tocolumn(y);
    end
    
    if size(y,1)==size(x,1) && size(y,2)==size(x,2)
        % pass
    elseif size(y,1)==size(x,1) && size(x,2)==1
%         x = repmat(x, [1,size(y,2)]);
    elseif size(y,2)==size(x,1) && size(y,1)==size(x,2)
        y = y';
    end
    
    %%% now x and y must have the same dimensions! %%%
    
    results = struct;
    
    for ii = 1:size(y,2)
        
        if size(x,2)==1
            xdata = x;
        else
            xdata = x(:,ii);
        end

        ydata = y(:,ii);
        
        if isempty(input.variances)
            vdata = nanvar(ydata).*ones(size(ydata,1),1);             
        elseif isscalar(input.variances)
            vdata = input.variances.*ones(size(ydata,1),1); 
        elseif size(input.variances,2)==1
            vdata = input.variances;
        elseif size(input.variances,1)==1 && size(input.variances,2)==size(y,2)
            vdata = ones(size(ytemp,1),1).*input.variances(1,ii);
        else
            vdata = input.variances(:,ii);
        end

        xdata = fillmissing(xdata, 'linear'); % there shouldn't be any NaNs in the xdata!
        
        bad_idx = isnan(ydata) | isnan(vdata);
                
        X = ones(length(xdata),1); 
        
        for k = 1:input.order
            X = [X xdata.^k];
        end
        
        for jj = 1:input.iterations

            ytemp = ydata(~bad_idx);
            vtemp = vdata(~bad_idx);
            Xtemp = X(~bad_idx,:); 
            
            C = 1./vtemp;
            XCX = Xtemp'*(C.*Xtemp); 
            
            results(ii).coeffs = XCX\(Xtemp'*(C.*ytemp));
            results(ii).model = print_model(input.order);
            
            warning('off', 'MATLAB:nearlySingularMatrix');
            warning('off', 'MATLAB:singularMatrix');

            y_model = X*results(ii).coeffs; % note we get the model for the entire dataset, including NaNs
            
            results(ii).coeffs = results(ii).coeffs';
            results(ii).variance = nanvar(ydata - y_model); 
            
            if isempty(input.variances) % vdata is calculated from the data, needs to be updated... 
                vdata = ones(size(ydata,1),1).*results(ii).variance; 
            end
            
            residuals = (ydata - y_model).^2./vdata;
            
            results(ii).chi2 = nansum(residuals); 
            
            results(ii).ndof = numel(ytemp)-input.order;
            
            % remove outliers 
            bad_idx = bad_idx + residuals>input.sigma.^2;

%             vtemp = vtemp./mean(vtemp,1,'omitnan').*results(ii).variance;
            
            if input.plotting
                
                if isempty(input.axes)
                    input.axes = gca;
                end
                
                plot(input.axes, xdata, ydata, '.', xdata, y_model, '-');
                
                pause(input.pause);
                
            end
            
            if nnz(bad_idx)==0 % further iterations don't change the values in xtemp/ytemp
                break;
            end
            
        end % for jj (iterations)
                
        results(ii).x = xdata;
        results(ii).y = ydata;
        results(ii).v = vdata;
        results(ii).ym = y_model;            
        results(ii).bad_idx = bad_idx; 
        
    end % for ii (fluxes)
    
end

function str = print_model(order)

    for ii = 1:order+1

        if ii==1
            str = 'y = coeff(1)';
        elseif ii==2
            str = [str ' + coeff(2)*x'];
        else
            str = sprintf('%s + coeff(%d)*x^%d', str, ii, ii-1);
        end

    end
    
end
    
    
    