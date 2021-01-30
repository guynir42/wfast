function I = gaussian(varargin)
% usage: gaussian(sigma_x,sigma_y=sigma_x,rot_frac=0,S=auto, norm=0)
% Generates an image of a 2D Gaussian at the center. 
% PARAMETERS:
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

    if nargin==0, help('util.shapes.gaussian'); return; end

    input = util.text.InputVars;
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
    [x,y] = util.shapes.make_grid(input.S, input.x_shift, input.y_shift, input.rot_frac); 
    
    I = exp(-0.5*((x./input.sigma_x).^2 + (y./input.sigma_y).^2));
    
    if input.norm==1
        I = I./util.stat.sum2(I);
    elseif input.norm==2
        I = I./sqrt(util.stat.sum2(I.^2));
    end
        
end