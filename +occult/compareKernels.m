function snr_frac = compareKernels(this_kernel, that_kernel, full_xcorr)
% Usage: snr_frac = compareKernels(this_kernel, that_kernel, full_xcorr)
% Calculate the fraction of S/N that remains when using the wrong kernel 
% (argument 2, "that_kernel") instead of filtering with the right kernel 
% (argument 1, "this_kerenl"). In both cases the S/N is calculated for an 
% input flux "this_kernel". 
%
% NOTE: "this_kernel" should be a column vector, "that_kernel" can be a 2D 
%       matrix, so each column is compared to "this_kernel", and the results
%       are returned as a row vector in "snr_frac".
%       The two inputs must have the same size in dimension 1. 
%
% The 3rd argument, "full_xcorr" (Default 1) calculates the full cross 
% correlation (with any shift). If set to 0, it will only calculate the 
% correlation when both kernels are perfectly centered (which is quicker). 
% 

    if nargin==0, help('occult.compareKernels'); return; end

    if nargin<3 || isempty(full_xcorr)
        full_xcorr = 1;
    end

    f1 = this_kernel - 1;
    f2 = that_kernel - 1;

    % now the normalized kernels
    k1 = f1./sqrt(sum(f1.^2));
    k2 = f2./sqrt(sum(f2.^2));

    if full_xcorr==0
        self_signal = sum(k1.*f1);
        cross_signal = sum(k2.*f1);
    else
        self_signal = max(filter2(f1,k1));
        cross_signal = max(filter2(f1,k2));
    end

    snr_frac = cross_signal./self_signal;

end