function f_new = resample(f, factor, varargin)
% Usage: f_new = resample(f, factor, varargin)

    if nargin==0, help('util.series.resample'); return; end
    
    S = size(f); 
    f = fillmissing(f, 'spline'); 
    
    F = fftshift(fft(f), 1); 
    
    F = util.img.fit2size(F, [ceil(S(1).*factor), S(2:end)]); 
    
    f_new = real(ifft(fftshift(F,1)))*factor; 
    
end