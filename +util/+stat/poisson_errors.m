function [lower, upper] = poisson_errors(lambda, p_value, precision_bytes)
% Usage: error_bounds = poisson_errors(lambda, p_value, precision_bytes=16)
% Calculate the upper and lower lambda values that have no more than
% p_value chance of generating the measured lambda value. 
% For example: if you measured lambda=10 counts, for all values below 
% 5.4254 there is less than 5% chance of getting 10 counts. 
%
% The input "precision_bytes" controls how many iterations the scan goes
% through before settling on a number. This roughly represents factors of 2
% in precision (relative to the input lambda!). 

    if nargin==0, help('util.stat.poisson_errors'); return; end
    
    if nargin<2 || isempty(p_value)
        p_value = 0.05; 
    end
    
    if nargin<3 || isempty(precision_bytes)
        precision_bytes = 16; 
    end
    
    if isempty(lambda)
        lower = [];
        upper = [];
    elseif ~isscalar(lambda)
        
        for ii = 1:numel(lambda)
            [lower(ii), upper(ii)] = util.stat.poisson_errors(lambda(ii), p_value, precision_bytes);
        end
        
        lower = reshape(lower, size(lambda));
        upper = reshape(upper, size(lambda));
        
    else

        % parse p_values for different input methods:
        if p_value < 0.5
            % do nothing
        elseif p_value>= 0.5 && p_value < 1 % we have been given confidence interval, instead of p-value
            p_value = 1 - p_value;
        elseif p_value>= 1 && p_value<50 % given p-value as percent points! 
            p_value = p_value/100;
        elseif p_value>=50 && p_value<100
            p_value = (100-p_value)/100;
        else
            error('Unreasonable p_value= %f. Input should be between 0 and 1 or between 1 and 100', p_value); 
        end

        % find lower bound
        low = 0;
        high = lambda;
        p = 1 - p_value; 

        for ii = 1:precision_bytes

            new_val = (high+low)/2; 

            if poissinv(p, new_val)>=lambda % too close to lambda
                high = new_val; % search the lower half
            else % far enough from lambda, we can get closer
                low = new_val; % search the upper half
            end

        end

        lower = low; 

        % find upper bound
        low = lambda;
        high = poissinv(1-p_value/4, lambda);
        p = p_value; 

        for ii = 1:precision_bytes

            new_val = (high+low)/2; 

            if poissinv(p, new_val)<=lambda % too close to lambda
                low = new_val; % search the upper half
            else % far enough from lambda, we can get closer
                high = new_val; % search the lower half
            end

        end

        upper = high; 

    end
    
end