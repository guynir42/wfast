function I = defocus_annulus(varargin)
% usage: I = defocus_annulus(r, minor_axis=r, sigma=2, dx=0, dy=0, rot_frac=0, S=auto, norm=0)
% Generates an image of an annulus apodized with Gaussian falloff, to look 
% similar to the defocused PSF when using a telescope with a central obscuration. 
%
% PARAMETERS:
%   -r: semi-major axis of the ellipse (or radius of the ring).
%   -minor_axis: if the ring is eccentric. Default is equal to r. 
%   -sigma: the gaussian falloff width parameter (default 2 pixels). 
%   -x_shift: move the center of the ring in the x direction.
%   -y_shift: move the center of the ring in the y direction.
%   -rot_frac: rotation angle in units of 0 to 1 (=90 degrees). 
%   -S: size of image (assumed square). Default is 2*r+5*sigma,
%       adjusted to be an odd-number.
%   -norm: choose normalization option. 
%          norm=0: not normalized (ring peak is 1);
%          norm=1: the sum of the shape is 1. 
%          norm=2: the sqrt of sum of squares of the shape is 1. 
%
% Parameters can be given in order or as keyword-value pairs...

    if nargin==0, help('util.shapes.ellipse'); return; end

    input = util.text.InputVars;
    input.input_var('r', [], 'a', 'radius', 'semi major axis');
    input.input_var('minor_axis', []); 
    input.input_var('sigma', 2, 'width'); 
    input.input_var('dx', 0, 'x_shift');
    input.input_var('dy', 0, 'y_shift');
    input.input_var('rot_frac', 0);
    input.input_var('S', [], 'size', 'imsize');
    input.input_var('norm', 0, 'normalization');
    input.use_ordered_numeric = 1;
    input.scan_vars(varargin{:});
    
    if isempty(input.r)
        error('Must supply a radius/semi-major axis');
    end
    
    if isempty(input.minor_axis)
        input.minor_axis = input.r;
    end
    
    if isempty(input.S)
        input.S = 2.*input.r + 5.*input.sigma;
        input.S = input.S+mod(input.S+1,2); % make sure it is an odd number
    end
    
    input.S = util.vec.imsize(input.S); % make sure it is a two element vector
    
    [x,y] = util.shapes.make_grid(input.S, input.dx, input.dy, input.rot_frac);
    
    % shorthands
    r = input.r;
    r2 = input.minor_axis;
    
    R = sqrt(x.^2./r.^2 + y.^2./r2.^2)*r;
    
    I = exp( -0.5*(R-r).^2./input.sigma.^2);
    
    if input.norm==1
        I = I./util.stat.sum2(I);
    elseif input.norm==2
        I = I./sqrt(util.stat.sum2(I.^2));
    end
    
end



