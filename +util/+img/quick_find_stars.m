function [table_props, I_reduced] = quick_find_stars(I, varargin)
% Usage: [table_props, I_reduced] = quick_find_stars(I, varargin)
% Finds all the stars in the image above some threshold, using regionprops
% 
% OPTIONAL ARGUMENTS:
% -number: Maximum number of stars to find. Default is Inf. 
% -dilate: How much to dilate each found star before subtracting it (default 9). 
% -psf_sigma: Estimate of the PSF gaussian width parameter (default 2).
% -mean: Estimate of the image mean (background). If empty will be calculated using im_stats.
% -std: Estimate of the image standard deviation. If empty will be calculated using im_stats.
% -sigma: Threshold for finding stars (after mean subtract and divide by std). Default 3.5.
% -saturation: Value of brightest pixel must be lower than this or else flag 1 is set. 
% -edges: Replace edges with NaN before finding stars. If empty will be dilate*2. 
%
% OUTPUTS: 
% -table_props: output from regionprops in a table, along with positions 
%               (centroids), flux and flag for each star. 
% I_reduced: the image with found stars replaced with NaNs. 
% 
% NOTE: flag values are given as:
% 1- saturated star
% 2- better fit to point source. 
% 3- better fit to extended source. 
% ...
%

    if nargin==0, help('util.img.quick_find_stars'); return; end

    input = util.text.InputVars;
    input.input_var('number', Inf, 'max_number', 'num_stars'); 
    input.input_var('dilate', 9);
    input.input_var('psf_sigma', 2);
    input.input_var('mean', []);
    input.input_var('std', []);
    input.input_var('sigma', 3.5);
    input.input_var('saturation', 5e6); % saturation by default for 100 images in a stack
    input.input_var('edges', []);
    input.scan_vars(varargin{:});
    
    if isempty(input.edges)
        input.edges = input.dilate.*2;
    end
    
    if isempty(input.mean) || isempty(input.std)
        
        [M,V] = util.img.im_stats(I);
        
        if isempty(input.mean)
            input.mean = M;
        end
        
        if isempty(input.std)
            input.std = sqrt(V);
        end
        
    end
    
    I = regionfill(I, isnan(I));
    I_reduced = I;
    
    if input.edges
        I_reduced = util.img.pad2size(util.img.crop2size(I_reduced, size(I)-input.edges.*2), size(I), NaN);
    end
    
    pos = [];
    flux = [];
    flag = [];
    table_props = table.empty;
    
    thresholds = input.sigma + [5 4 3 2 1 0];
    
    k = util.img.gaussian2(input.psf_sigma, 'norm', 2);
    k_half = util.img.gaussian2(input.psf_sigma./2, 'norm', 2, 'size', size(k));
    k_twice = util.img.gaussian2(input.psf_sigma.*2, 'norm', 2, 'size', size(k));

    for ii = 1:length(thresholds)
        
        I2 = I_reduced;
        I2(isnan(I2)) = 0;
        
        If = filter2(k, I_reduced)./sqrt(util.stat.sum2(k)); 
        If_half = filter2(k_half, I_reduced)./sqrt(util.stat.sum2(k_half));
        If_twice = filter2(k_twice, I_reduced)./sqrt(util.stat.sum2(k_twice));

        BW = (If-input.mean)./input.std>=thresholds(ii);

        BW_dilated = imdilate(BW, ones(input.dilate));

        T = regionprops('table', BW, I, 'WeightedCentroid', 'PixelValues', 'PixelIdxList'); 

        N = height(T);

        if N>0
            
            F = cellfun(@sum, T.PixelValues);

            P = T.WeightedCentroid;

            MX = zeros(N,1);
            MX_half = zeros(N,1);
            MX_twice = zeros(N,1);

            for jj = 1:N
                MX = max(If(T.PixelIdxList{jj}));
                MX_half = max(If_half(T.PixelIdxList{jj}));
                MX_twice = max(If_twice(T.PixelIdxList{jj}));
            end

            f = zeros(size(P,1),1);    
            f(MX<MX_half) = 2; % point source 
            f(MX<MX_twice) = 3; % extended object
            f(MX>input.saturation) = 1; % saturated star 

            pos = vertcat(pos, P);
            flux = vertcat(flux, F);
            flag = vertcat(flag, f);
            table_props = vertcat(table_props, T);

            if height(table_props)>input.number
                break; % shorten the run if we already found enough stars...
            end
            
            I_reduced(BW_dilated) = NaN;
        
        end
        
    end
    
    T2 = table(flux, flag);
    
    table_props = horzcat(T2, table_props);
    
    table_props = sortrows(table_props, 1, 'descend');
    
    if ~isempty(table_props)

        table_props.Properties.VariableNames{'WeightedCentroid'} = 'pos';
        table_props = [table_props(:,'flux') table_props(:,'flag') table_props(:,'pos') table_props(:, 'PixelValues'), table_props(:, 'PixelIdxList')];

        if height(table_props)>input.number
            table_props = table_props(1:input.number, :);
        end

    end
    
end