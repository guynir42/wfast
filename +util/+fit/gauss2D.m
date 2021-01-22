function results = gauss2D(Im, varargin)
% Usage: results = gauss2D(Im, varargin)
% Fit a 2D gaussian to an image. 
% Use a first-guess linear fit to the log of the image, then use non-linear
% minimization to find the best fit with more advanced models. 
% 
% Input: Im is the target image. Should be a small cutout with one source. 
%        If multiple images are given in a 3D (or higher) matrix, the same
%        algorithm is run in a loop over each image. 
% 
% Output: results is a structure with some model parameters, the resulting 
%         model image, the difference image, the Chi2 score, etc. 
%         If the input is multi-image, the function returns a struct array
%         where the dimensions agree with the 3rd and higher dims of Im. 
%
% OPTIONAL ARGUMENTS:
%   -
%
% 

    if nargin==0, help('util.fit.gauss2D'); return; end
    
    input = util.text.InputVars;
    input.input_var('linear', true); % use a pre-fit with a simple gaussian model and linear fit to log(Im)
    input.input_var('search', true); % use a non-linear fit with more complicated model
    input.input_var('threshold', 1); % pixel value below which we replace the values with NaN
    input.input_var('aperture', [], 'radius'); % aperture radius (in pixels) beyond which we place NaNs (default is not to use). 
    input.input_var('generalized_gaussian', false); % use a generalized gaussian instead of a regular one
    input.input_var('plot', false, 'use_plot'); % plot the fitting process
    input.input_var('pause', 0, 'pause_duration'); % pause between fitting iterations
    input.input_var('axes', [], 'axis'); % axes to draw into, when plotting (default is gca())
    input.scan_vars(varargin{:}); 
    
    Im_all = Im; 
    S_all = size(Im_all);
    if ndims(Im)>2
        N = prod(S_all(3:end)); % number of images total
    else
        N = 1;
    end
    
    for ii = 1:N
    
        Im = Im_all(:,:,ii); % loop over each image separately
        
        Im(Im<input.threshold) = NaN;
        if ~isempty(input.aperture)
            circ = util.img.ellipse(input.aperture, 'size', size(Im)); 
            Im(~circ) = NaN; 
        end
        
        if input.linear
            [coeffs, model_im] = linearFitGaussian(Im); % quick linear fitter
        else
            coeffs.amp = 1; 
            coeffs.x_offset = 0;
            coeffs.y_offset = 0;
            coeffs.sig_x = 0; 
            coeffs.sig_y = 0;
            coeffs.rot = 0;

            model_im = [];
        end
        
        results(ii).coeffs = coeffs;
        results(ii).model_image = model_im;
        results(ii).diff_image = Im - model_im; 
    
    end
    
    if ndims(Im)>2 % if not, results is a scalar struct
        results = reshape(results, S(3:end)); % make sure the dimensions follow the image matrix dimensions
    end
    
end

function [coeffs, model_im] = linearFitGaussian(Im)
% Some algebra: we want to fit Im to a gaussian in 2D:
%  f(x,y | amp,mu_x,mu_y,L_x,L_y,L_xy) = amp*exp[-0.5*(L_xx*(x-mu_x)^2 ...
%     + L_yy*(y-mu_y)^2 + 2*L_xy*(x-mu_x)(y-mu_x)
%
% Taking the log of this function we get:
%  log(f) = log(amp) - 0.5*[L_xx*(x^2 - x*2*mu_x + mu_x^2)...
%     + L_yy*(y^2 - y*2*mu_y + mu_y^2) ...
%     +2*L_xy*(x*y + x*mu_y -y*mu_x + mu_x*mu_y]
% 
% We can parameterize using:
%  f(x,y) = C1 + C2*x + C3*y + C4*x*y + C5*x^2 + C6*y^2 
% with these parameters: 
%  C1 = log(amp) - 0.5*L_xx*mu_x^2 - 0.5*L_yy*mu_y^2 - L_xy*mu_x*mu_y
%  C2 = L_xx*mu_x + L_xy*mu_y
%  C3 = L_yy*mu_y + L_xy*mu_x
%  C4 = -L_xy
%  C5 = -0.5*L_xx
%  C6 = -0.5*L_yy

    [X,Y] = meshgrid((1:size(Im,2))-floor(size(Im,2)/2)-1, (1:size(Im,1))-floor(size(Im,1)/2)-1); 
    
    IL = log(Im); 
    idx = ~isnan(IL); % good indices from this matrix
    measurements = IL(idx); 
    weights = (Im(idx)).^2; % set the weights to the flux values? 
    
    X = X(idx); % only keep X values of good pixels
    Y = Y(idx); % only keep Y values of good pixels
    
    A = [ones(length(measurements),1), X, Y, X.*Y, X.^2, Y.^2]; % design matrix
    
    C = lscov(A, measurements, weights); % solve for the coefficients C
    
    % solve for physical parameters:
    L_xx = -2*C(5);
    L_yy = -2*C(6);
    L_xy = -C(4);    
    mu_y = (C(3)/L_yy - C(2)*L_xy/L_xx/L_yy)/(1-L_xy^2/L_xx/L_yy); 
    mu_x = (C(2) - L_xy*mu_y)/L_xx;     
    amp = exp( C(1) + 0.5*L_xx*mu_x^2 + 0.5*L_yy*mu_y^2 + L_xy*mu_x*mu_y );
    
    % find the rotation and major/minor axes of L (the inverse covariance)
    L = [L_xx, L_xy; L_xy, L_yy]; 
    [R, E] = eig(L); % get the rotation matrix and eigenvector matrix for L
    
    coeffs.amplitude = amp;
    coeffs.x_offset = mu_x;
    coeffs.y_offset = mu_y;
    coeffs.sigma_x = 1./sqrt(E(4)); % inverse covariance eigenvalue
    coeffs.sigma_y = 1./sqrt(E(1)); % inverse covariance eigenvalue
    coeffs.angle = 180 + atan2d(R(4), R(2)); 
    
    % now recreate the gaussian based on the parameters found
    model_im = amp.*util.img.gaussian2('sigma_x', coeffs.sigma_x, 'sigma_y', coeffs.sigma_y, ...
        'x_shift', coeffs.x_offset, 'y_shift', coeffs.y_offset, 'rot_frac', coeffs.angle/90, 'size', size(Im)); 
    
end

function sum_of_differences = helper_function_regular(beta, Im, use_plot, pause_duration, axes)
    
end

