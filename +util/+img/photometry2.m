% Usage: [outputs, arrays] = photometry2(cutouts, varargin)
% Calculate all sorts of photometric measurements on image cutouts. 
%
% By default, this calculates raw photometry (sum of the entire square), 
% aperture photometry on a few radii, forced photometry (on the largest
% radius using the average position of all stars), and gaussian PSF
% photometry, assuming some gaussian width. 
% 
% The position of the offsets of each star in its cutout is calculated by
% first taking the centroids of the background subtracted sum on the whole
% square (raw-photometry). 
% If use_centering=1 then an aperture is moved onto the offset position and
% used to produce new centroids. 
% If use_gaussian=1 then PSF photometry is performed iteratively, each time
% improving the centroid positions. 
% Those positions are used to place a few concentric apertures, calculating
% the photometry in that position for each radius. 
% Finally, the flux-averaged offsets are calculated from the best position
% estimate (which is chosen to be gaussian, aperture or raw, in decending
% order of accuracy), and used to do forced photometry on all cutouts with
% the largest aperture placed at the average offsets. 
% In all cases the background is estimated from an annulus centered around 
% at the position of the PSF/aperture. For raw-photometry it is centered
% around the middle of the cutout. 
%
% Output: Each type of photometry products are returned in a separate
%         struct, where each field contains a matrix of results for each
%         star and each cutout. 
%         The different products are:
%           *flux: not background subtracted! 
%           *area: how many pixels in the aperture/PSF, removing bad pixels
%           *error: estimate of the photometric error per sample. 
%           *background: per pixel, to subtract do flux-area*background. 
%           *variance: per pixel, measured from the background annulus. 
%           *offset_x/y: in pixels, relative to the center of the cutout. 
%           *width: average 2nd moment of the PSF, equivalent to gaussian
%                   "sigma" width parameter. 
%           *bad_pixels: that are inside the effective aperture.
%                        (what does this mean for gaussians?)
%           *flag: set to 1 if there is a problem with the centroids 
%                  (e.g., negative width, very large offsets). 
%
% Optional Arguments:
%       *
%
%
%
%
