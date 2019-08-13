function results = polyfit(x,y,varargin)
% Usage: results = polyfit(x,y,varagin)

    if nargin==0, help('util.fit.polyfit'); return; end

    input = util.text.InputVars;
    input.input_var('order', 2); % polynomial order
    input.input_var('sigma', 3); % above this threshold is considered outliers
    input.input_var('iterations', 3); % how many iterations of removing outliers
    input.input_var('variances', []); % assume no covariance, only different variance per sample
    input.input_var('double', 0); % use double precision
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
            xtemp = x;
        else
            xtemp = x(:,ii);
        end

        ytemp = y(:,ii);
        
        if isempty(input.variances)
            vtemp = var(ytemp,[],1,'omitnan').*ones(size(ytemp,1),1); 
        elseif isscalar(input.variances)
            vtemp = input.variances.*ones(size(ytemp,1),1); 
        elseif size(input.variances,2)==1
            vtemp = input.variances;
        elseif size(input.variances,1)==1 && size(input.variances,2)==size(y,2)
            vtemp = ones(size(ytemp,1),1).*input.variances(1,ii);
        else
            vtemp = input.variances(:,ii);
        end

        nan_idx = isnan(xtemp) | isnan(ytemp) | isnan(vtemp);

        xtemp(nan_idx) = [];
        ytemp(nan_idx) = [];
        vtemp(nan_idx) = [];
        
        for jj = 1:input.iterations
            
            % reasses the variance 
            
            X = ones(length(xtemp),1); 
            for k = 1:input.order
                X = [X xtemp.^k];
            end
            
            C = 1./vtemp;
            XCX = X'*(C.*X); 
            
            results(ii).coeffs = XCX\(X'*(C.*ytemp));
            results(ii).model = print_model(input.order);
%             results(ii).coeffs = inv(XCX)*(X'*(C.*ytemp));
            
            warning('off', 'MATLAB:nearlySingularMatrix');
            warning('off', 'MATLAB:singularMatrix');

            y_model = X*results(ii).coeffs; 

            residuals = (ytemp - y_model).^2./vtemp;
            
            results(ii).chi2 = sum(residuals); 
            
            results(ii).ndof = numel(ytemp)-input.order;
            
            results(ii).variance = var(ytemp - y_model); 

            % remove outliers 
            bad_idx = residuals>input.sigma.^2;

            xtemp(bad_idx) = [];
            ytemp(bad_idx) = [];
            vtemp(bad_idx) = [];
            y_model(bad_idx) = [];

            results(ii).x = xtemp;
            results(ii).y = ytemp;
            results(ii).v = vtemp;
            results(ii).ym = y_model;            
            
            vtemp = vtemp./mean(vtemp,1,'omitnan').*results(ii).variance;
            
            if input.plotting
                
                if isempty(input.axes)
                    input.axes = gca;
                end
                
                plot(input.axes, xtemp, ytemp, '.', xtemp, y_model, '-');
                
                pause(input.pause);
                
            end
            
        end
        
    end
    
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
    
    
    