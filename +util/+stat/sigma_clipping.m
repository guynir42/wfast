function [mu, sigma, values, N_exc] = sigma_clipping(values, varargin)
% Usage: [mu, sigma, values, N_exc] = sigma_clipping(values, varargin)
% Calculate a fit to the distribution of "values", iteratively removing
% outliers and refitting. 
%
% Input: a vector of values to fit to a distribution. 
% Outputs: -the mean of the final, fitted distribution. 
%          -the standard devaition of the final distribution. 
%          -a vector of values after removing outliers. 
%          -the number of excluded values. 
%
% OPTIONAL ARGUMENTS:
%   *sigma: how many standard deviations from the mean we use to clip values. 
%           Default is 5. 
%   *iterations: how many maximal iterations should we do. The function
%                quits if no new values are excluded, or when reaching this
%                parameter value. Default is 5. 
%   *distribution: what to fit to the histogram of values. Default is
%                  "gauss" but can also handle "weibull" (extreme value). 
%   *plot: display the fitting process for debugging (default false). 
%   *axes: what axes to plot to. Default is "gca". 
%   *pause: add additional time between iteration when plotting. Default 0.
%   

    import util.text.*;
    import util.stat.*;
    
    if nargin==0, help('util.stat.sigma_clipping'); return; end
    
    % add input checks on "values"
    values = values(:);
    values(isnan(values)) = [];
    
    Nsigma = 5;
    Niter_max = 5;
    dist = 'gauss';
    use_plot = 0;
    ax = [];
    pause_length = 0;
    
    for ii = 1:2:length(varargin)
        key = varargin{ii};
        val = varargin{ii+1};
        
        if cs(key, {'sigma', 'nsigma', 'num_sigma', 'number sigma'})
            Nsigma = val;
        elseif cs(key, 'iterations')
            Niter_max = val;
        elseif cs(key, 'distribution')
            dist = val;
        elseif cs(key, 'plotting')
            use_plot = parse_bool(val);
        elseif cs(key, {'axis', 'axes'})
            ax = val;
        elseif cs(key, 'pause length')
            pause_length = val;
        end
        
    end
    
    if isempty(ax) && use_plot
        ax = gca;
    end
    
    N_prev = numel(values);
    values = values(~isnan(values));
    
    % first iteration:
    mu = median(values, 'omitnan');
%     sigma = fwhm(values);
    sigma = std(values, 'omitnan');
%     N = sum(values>mu-Nsigma*sigma & values<mu+Nsigma*sigma);
    
    values = values(values>mu-Nsigma*sigma & values<mu+Nsigma*sigma);
    N = numel(values);
    N_exc = N_prev - N; % number of excluded values
    
    if isnan(mu)
        return; % silently give back a NaN if the inputs are all NaN!
    end
    
    if N < 4
        return;
    end
    
    for ii = 1:Niter_max
        
        mu_prev = mu;
        sigma_prev = sigma;
        N_prev = N;
        
        if cs(dist, 'gauss')
%             phat = mle(values, 'logpdf', @(x,m,s) -0.5*log(2*pi*s.^2) -0.5*(x-m).^2/s.^2, 'start', [mu_prev,sigma_prev],...
%                 'LowerBound',[mu_prev-Nsigma*sigma_prev, 0], 'UpperBound', [mu_prev+Nsigma*sigma_prev, Nsigma*sigma_prev]);
            phat = fitdist(values, 'normal');
            mu = phat.mu;
            sigma = phat.sigma;
        elseif cs(dist, {'weibull', 'max values'})
            phat = evfit(-double(values)); 
            
            mu = -phat(1);
            sigma = phat(2);
        end
        
        
        if use_plot
            disp(['plotting iteration ' num2str(ii)]);
            step = (max(values)-min(values))/sqrt(N);
            if all(values==round(values)) && step<1
                step = 1;
            end
                
            x = min(values):step:max(values);
            x2 = min(values):step/4:max(values);
            histogram(ax, values, x);
            
            if cs(dist, 'gauss')
                line(ax, x2, step*N*1/sqrt(2*pi*sigma.^2)*exp(-0.5*((x2-mu)/sigma).^2), 'color', 'red');
            elseif cs(dist, {'weibull', 'max values'})
                line(ax, x2, step*N*evpdf(-x2, -mu, sigma), 'color', 'red');
            end
            
            pause(pause_length);
            
        end
        
        keep_idx = values>mu-Nsigma*sigma & values<mu+Nsigma*sigma;
        
        values = values(keep_idx);
                
        N = numel(values);        
        
        N_exc = N_exc + N_prev - N;

        if N==N_prev || N<4
            return;
        end
        
    end
        
end