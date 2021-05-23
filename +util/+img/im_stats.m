function [M,V] = im_stats(Im, varargin)
% Usage: [M,V] = im_stats(Im, varargin)
% Gets basic statistics on a 2D image by cutting it into smaller tiles
% (using "jigsaw" function) and calculating the average median and variance
% values. This is a quick and dirty way to estimate the background and
% noise in an image even if it has stars in it. 
% 
% OPTIONAL ARGUMENTS:
%   *method: how to calculate the mean/variance in each tile. 
%            Choose median (default) for simple median2/var2 calculation. 
%            Choose clipping to use sigma clipping on each tile. 
%   
%   *output: choose if you want to get a single median value for the whole
%            image ("median" is the default) or an interpolated map of values
%            for each point in the image ("map"). 
%
% The rest of the arguments given are passed to "jigsaw" as they are. 
% This includes:, tile size, padding, overlap, squeeze, partial. 
% 
    
    import util.text.cs;

    if nargin==0, help('util.img.im_stats'); return; end

    input = util.text.InputVars;
    input.input_var('method', 'median'); 
    input.input_var('output', 'median'); 
    input.scan_vars(varargin{:}); 
    
    if isempty(Im)
        error('Got an empty input!'); 
    end
    
    [C, P, D] = util.img.jigsaw(Im, 'partial', 0, varargin{:}); % get cutouts, positions of centers, and dimensions
    
    if cs(input.method, 'median')
        
        bias_values = util.stat.median2(C);
%         var_values = util.stat.var2(C);
        var_values = util.stat.rstd(reshape(C, [size(C,1).*size(C,2), 1, size(C,3)])).^2; 
    
    elseif cs(input.method, 'clipping', 'sigma_clipping')
        
        for ii = 1:size(C,3)
            [bias_values(ii), var_values(ii)] = util.stat.sigma_clipping(C(:,:,ii), 'sigma', 2.5, 'iter', 3); 
            var_values(ii) = var_values(ii).^2;
        end
        
    else
        error('Unknown "method" option: "%s". Use "median" or "clipping" instead...', input.method); 
    end
    
    if cs(input.output, 'median')
    
        M = median(bias_values, 'omitnan');
        V = median(var_values, 'omitnan');

    elseif cs(input.output, 'map')
        
        bias_values = reshape(bias_values, D)'; % reshape to a 2D grid of values for each cutout center point
        var_values = reshape(var_values, D)'; % reshape the variance also
        
        [X,Y] = meshgrid(unique(P(:,1)), unique(P(:,2))); % the X/Y points of the cutout centers
        
        [xq,yq] = meshgrid(1:size(Im,2), 1:size(Im,1)); % grid of query points is for each point in the original image
        
        M = interp2(X, Y, bias_values, xq, yq, 'spline'); % get an interpolation of all grid points
        V = interp2(X, Y, var_values, xq, yq, 'spline'); % get an interpolation of all grid points
        
    else
        error('Unknown "output" option: "%s". Choose "median" or "map" instead...', input.output); 
    end
    
end