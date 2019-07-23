function values_conv = convolution(kernels, values, varargin)
% Usage: values_conv = convolution(kernels, values, varargin)
% Run convolution or cross-correlation of "kernels" with "values". 
% Both inputs can be any size vector or matrix, but any non-singleton dimensions
% beyond dim1 must match (i.e., they are both padded to the same length in 
% dimension 1 and then multiplied. This also gives the dimensions of the output. 
%
% This function runs in Fourier space and requires a lot of memory if these 
% are big matrices, especially if one has large dim2 and the other dim3 
% (the output and intermediate matrices will have dim2*dim3 size). 
%
% OPTIONAL ARGUMENTS:
%  *cross: cross-correlation instead of convolution (i.e., use complex conjugate in Fourier space). Default is false. 
%  *crop: what is the final region of the result. 
%         "same" (default) outputs the same size, "valid" gives only the 
%         region the has full overlap, "full" gives the joint size with mutual padding.
%  *debug_bit: Level of verbosity (default 0). 
%  *mem_max: Maximum available memory (in GBs) for multiplying big matrices. 
%            If the calculation requires more than this, it is done in a 
%            very slow loop... 

    input = util.text.InputVars;
    input.input_var('cross', 0, 'conj'); 
    input.input_var('crop', 'same'); 
    input.input_var('debug_bit', 0); 
    input.input_var('mem_max', 10, 'memory max', 'memory limit'); 
    input.scan_vars(varargin{:});
    
    L = size(kernels,1) + size(values,1) - 1; % length of dim1 for both kernels and values
            
    Sk = size(kernels);
    
    k = util.img.pad2size(kernels, [L Sk(2:end)]);
    
    kernels_fft = fft(k);
    
    if input.cross
        kernels_fft = conj(kernels_fft);
    end

    Sf = size(values);

    % if values ends up having many dimensions (above 4) we will need to patch the pad2size function 
    f = util.img.pad2size(values, [L Sf(2:end)]); % keep all dimensions of values except the first, which is padded
    
    values_fft = fft(f);

    M = L*size(kernels_fft,2)*size(values_fft,3)*4/1024^3;
    
    if input.debug_bit>1, disp(['memory usage: ' num2str(M) ' GBs']); end

    if M<input.mem_max
        values_conv = real(fftshift(ifft(kernels_fft.*values_fft),1));
    else
        error('Memory required is too large: %f GBs. Use "mem_max" to set it higher...', M);
%         % maybe add preallocation? 
%         for ii = 1:size(values_fft,3)
%             values_conv(:,:,ii) = real(fftshift(ifft(kernels_fft.*values_fft(:,:,ii)),1));
%         end
    end

    if util.text.cs(input.crop, 'same')
        values_conv = util.img.crop2size(values_conv, [Sf(1), Sk(2), Sf(3:end)]); % crop back to original dimensions...
    elseif util.text.cs(input.crop, 'valid')
        
        S = size(values_conv);
        S(1) = length(values)-length(kernels)+1;
        values_conv = util.img.crop2size(values_conv, S);
%         values_conv = values_conv(1:S(1), :, :,:);
    elseif util.text.cs(input.crop, 'full')
        % pass
    else
        error('Unknown "crop" value "%s". Try "same" or "valid" or "full".', input.crop);
    end


end