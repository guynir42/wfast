function flux_det = detrend(flux, varargin)
% Usage: flux_det = detrend(flux, varargin)
% Remove trends from fluxes. 
% 
% INPUT: flux values (lightcurves). If a single vector, it can be row or 
%        column. If more dimensions exist it will treat dim1 as time.
%
% OPTIONAL ARGUMENTS:
%   -order: the polynomial order used in fitting the trend. Default 1. 
%
%
% OUTPUT: The output is the same size as the input. 

    if nargin==0, help('util.series.detrend'); return; end
    
    input = util.text.InputVars; 
    input.input_var('order', 1); % do we need this? 
    input.scan_vars(varargin{:}); 
    
    S_original = size(flux); 
    
    if isvector(flux)
        flux = util.vec.tocolumn(flux);
    else
        flux = flux(:,:); % linearize higher dimensions
    end
    
    fr = util.fit.polyfit(1:size(flux,1), flux, 'double', 1, 'order', 1, 'quick', 1, varargin{:}); 
    
    flux_det = flux - [fr.ym]; 
    
    flux_det = reshape(flux_det, S_original); 
    
end