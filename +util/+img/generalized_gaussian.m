function I = generalized_gaussian(varargin)
% Usage: I = generalized_gaussian(varargin)
% A more generic functional shape that better represents a defocused PSF. 
%
% Formula ( p^(1-1/p) ) / (2 sigma Gamma(1/p) ) * exp (-1/p * |x-mu|^p / sigma^p)
% is taken from Claxton & Staunton 2008
%    
% PARAMETERS:
%   -power: the "hardness" parameter of the generalized gaussian. 
%           When p=2 we get a regular gaussian. Increasing p will make the 
%           gaussian sharper, with p->inf giving a "pillbox" shape. 
%   -sigma_x: the width sigma parameter in the X direction. 
%   -sigma_y: the width sigma in Y. Default is same as X. 
%   -rot_frac: rotation angle in units of 0 to 1 (=90 degrees). 
%   -x_shift: move the center of the Gaussian in the x direction;
%   -y_shift: move the center of the Gaussian in the y direction;
%   -S: size of image (assumed square). Default is ceil(max(sigma_x,sigma_y)*10),
%       adjusted to be an odd-number.
%   -norm: choose normalization option. 
%          norm=0: not normalized (peak is 1);
%          norm=1: the sum of the Gaussian is 1. 
%          norm=2: the sqrt of sum of squares of Gaussian is 1. 
%
% Parameters can be given in order or as keyword-value pairs...

    if nargin==0, help('util.img.gaussian2'); return; end

    input = util.text.InputVars;
    input.input_var('power', 2, 'hardness'); 
    input.input_var('sigma_x', [], 'CX');
    input.input_var('sigma_y', [], 'CY');
    input.input_var('x_shift', [], 'dx');
    input.input_var('y_shift', [], 'dy');
    input.input_var('rot_frac', 0);
    input.input_var('S', [], 'size', 'imsize');
    input.input_var('norm', 0, 'normalization');
    input.use_ordered_numeric = 1;
    input.scan_vars(varargin{:});
    
    if isempty(input.sigma_x) && isempty(input.sigma_y)
        error('Must supply a width of the Gaussian...');
    elseif ~isempty(input.sigma_x) && isempty(input.sigma_y)
        input.sigma_y = input.sigma_x;
    elseif isempty(input.sigma_x) && ~isempty(input.sigma_y)
        input.sigma_x = input.sigma_y;
    end
        
    if isempty(input.S)
        input.S = ceil((max(input.sigma_x, input.sigma_y)*10));
        input.S = input.S+mod(input.S+1,2);
    end
    
    input.S = util.vec.imsize(input.S); 
    
    [x,y] = meshgrid(-floor((input.S(2))/2):floor((input.S(2)-1)/2), -floor((input.S(1))/2):floor((input.S(1)-1)/2));
    
    if isempty(input.rot_frac)
        x2 = x;
        y2 = y;
    else
        x2 = +x*cos(pi/2*input.rot_frac)+y*sin(pi/2*input.rot_frac);
        y2 = -x*sin(pi/2*input.rot_frac)+y*cos(pi/2*input.rot_frac);
    end
    
    if ~isempty(input.x_shift)
        x2 = x2 - input.x_shift;
    end
    
    if ~isempty(input.y_shift)
        y2 = y2 - input.y_shift;
    end
    
%     I = exp(-0.5*((x2./input.sigma_x).^2 + (y2./input.sigma_y).^2)); % simple gaussian... 
    
%   ( p^(1-1/p) ) / (2 sigma Gamma(1/p) ) * exp (-1/p * |x-mu|^p / sigma^p)
%     I = input.power^(1-1/input.power)./(2*input.sigma.*gamma(1/input.power)).*exp(-1/input.power*(abs(x2./input.sigma_x).^input.power + abs(y2./input.sigma_y).^input.power));
    I = exp(-1/input.power*(abs( (x2./input.sigma_x).^2 + (y2./input.sigma_y).^2 ).^(input.power/2))); % no need for normalizations! 
    
    if input.norm==1
        I = I./util.stat.sum2(I);
    elseif input.norm==2
        I = I./sqrt(util.stat.sum2(I.^2));
    end
        
end


% Citation: 
% @article{Claxton:08,
% author = {Christopher D. Claxton and Richard C. Staunton},
% journal = {J. Opt. Soc. Am. A},
% keywords = {Image detection systems; Optical transfer functions; Range finding; Propagating methods ; Testing ; Optical engineering; Edge detection; Image processing; Imaging noise; Low light levels; Machine vision; Spatial frequency},
% number = {1},
% pages = {159--170},
% publisher = {OSA},
% title = {Measurement of the point-spread function of a noisy imaging system},
% volume = {25},
% month = {Jan},
% year = {2008},
% url = {http://josaa.osa.org/abstract.cfm?URI=josaa-25-1-159},
% doi = {10.1364/JOSAA.25.000159},
% abstract = {The averaged point-spread function (PSF) estimation of an image acquisition system is important for many computer vision applications, including edge detection and depth from defocus. The paper compares several mathematical models of the PSF and presents an improved measurement technique that enables subpixel estimation of 2D functions. New methods for noise suppression and uneven illumination modeling were incorporated. The PSF was computed from an ensemble of edge-spread function measurements. The generalized Gaussian was shown to be an 8 times better fit to the estimated PSF than the Gaussian and a 14 times better fit than the pillbox model.},
% }