classdef MicroFlare < handle

    properties(Transient=true)
        gui;
    end
    
    properties % inputs/outputs
        
        filename; % where the flare was detected
        folder; % what folder the flare was in
        serial; % index of the total number of flares detected in this run
        file_index; % index in the run
        frame_index; % what frame did the peak appear in
        cut_index; % which cutout (star) from the file cutouts is this flare coming from
        pos; % x and y position of the peak pixel in the full-frame image
        peak; % value of the peak pixel
        pixel_var; % variance of underlying pixel (from calibration file)
        cutouts; % cutout images around the peak and all frames before/after it in that batch
        
        % from photometry2 or from rough estimates on the whole cutout
        timestamps;
        flux;
        error;
        background;
        area;
        offset;        
        width;
        bad_pixels;
        
        mean; % mean flux calculated outside the peak region
        std; % RMS flux calculated outside the peak region
        
        back_mean; % mean of images
        back_std; % RMS of images
        
        image; % single cutout at the frame of the peak
        num_peaks; % how many distinct peaks can we find in the one image
        num_frames; % how many continuous frames we have around the peak
        num_pixels; % how many pixels above the threshold do we have 
        
        fwhm; % calculted from the image, using dedicated FWHM function (e.g., using filters)
        
        type = ''; % can be cosmic ray, flare, satellite, pixel... 
        
    end
    
    properties % switches/controls
        
        pixel_thresh = 256; % threshold of the brightest pixel, used for triggering 
        region_thresh = 5; % each pixel is counted if it surpasses this many multiples of  image noise 
        frac_thresh = 0.1; % threshold of pixels around the brightest pixel, as fraction of peak, for finding number of pixels
        flux_thresh = 7.5; % threshold for peak intensity (of the flux, not of single pixel)
        side_thresh = 3; % threshold for flux in frames adjacent to the peak
        interval_peak = 2; % how many frames in either direction from peak should be excluded for calculating the mean/std of the flux
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        
        
    end
    
    properties(Hidden=true)
       
        version = 1.00;
        
    end
    
    methods % constructor
        
        function obj = MicroFlare(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'img.MicroFlare')
                if obj.debug_bit>1, fprintf('MicroFlare copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit>1, fprintf('MicroFlare constructor v%4.2f\n', obj.version); end
            
            end
            
        end
        
    end
    
    methods % reset/clear
        
    end
    
    methods % getters
        
    end
    
    methods % setters
        
    end
    
    methods(Static=true) % utilities
        
        function obj = readStruct( st)
            
            if ~isstruct(st)
                error('Must supply a struct array to this function! Instead got a "%s" type object.', class(st)); 
            end
            
            list = properties(img.MicroFlare); 
            
            for ii = 1:length(st)
                obj(ii) = img.MicroFlare; 
                for jj = 1:length(list)
                    if isfield(st(ii), list{jj})
                        obj(ii).(list{jj}) = st(ii).(list{jj}); 
                    end
                end
            end
            
        end
        
    end
    
    methods % calculations
        
        function calculate(obj, varargin)
            
            if numel(obj)>1
                
                for ii = 1:numel(obj)
                    obj(ii).calculate(varargin{:}); % recursively call this for individual objects in the array
                end
                
                return;
                
            end
            
            input = util.text.InputVars;
            input.input_var('pixel_thresh', [], 'pixel_threshold'); 
            input.input_var('area_thresh', [], 'area_threshold'); 
            input.input_var('flux_thresh', [], 'flux_threshold'); 
            input.input_var('side_thresh', [], 'side_threshold'); 
            input.input_var('interval_peak', []); 
            input.input_var('debug_bit', []); 
            input.scan_vars(varargin{:}); 
            
            for ii = 1:length(input.list_added_properties)
                if ~isempty(input.(input.list_added_properties{ii})) && isprop(obj, input.list_added_properties{ii})
                    obj.(input.list_added_properties{ii}) = input.(input.list_added_properties{ii}); 
                end
            end
            
            obj.image = obj.cutouts(:,:,obj.frame_index); 
            
%             indices = obj.frame_index;
%             for ii = 1:obj.interval_peak
%                 
%                 idx = obj.frame_index - ii;
%                 if idx<1, break; end
%                 
%                 indices = [indices idx]; 
%                     
%             end
%             
%             for ii = 1:obj.interval_peak
%                 
%                 idx = obj.frame_index + ii;
%                 if idx>length(obj.flux)
%                     break;
%                 end
%                 
%                 indices = [indices idx]; 
%                     
%             end
%             
%             indices = sort(indices);
%             flux_bg = obj.flux;
%             flux_bg(indices) = NaN;
%             
%             obj.mean = nanmean(flux_bg);
%             obj.std = nanstd(flux_bg); 
%             
%             flux_norm = (obj.flux-obj.mean)./obj.std;
            
            flux_norm = obj.getNormFlux; 

            % how many frames above the threshold
            counter = 1;
            for ii = 1:length(flux_norm)
                idx = obj.frame_index - ii;
                
                if idx<1, break; end
                if flux_norm(idx)<obj.side_thresh, break; end
                
                counter = counter + 1;
                
            end
            
            for ii = 1:length(flux_norm)
                idx = obj.frame_index + ii;
                
                if idx>length(flux_norm), break; end
                if flux_norm(idx)<obj.side_thresh, break; end
                
                counter = counter + 1;
                
            end
            
            obj.num_frames = counter;
            
            obj.back_mean = nanmedian(util.stat.mean2(obj.cutouts)); 
            obj.back_std = nanmedian(util.stat.std2(obj.cutouts)); 
            
            % how many pixels in the image are above the threshold
            BW = (obj.image-obj.back_mean)./obj.back_std>obj.region_thresh; % black/white image of pixels above the threshold
            N = nnz(BW); 
            
            if N<=1
                obj.num_pixels = 1;
            else
                
                C = bwconncomp(BW); 
                [X,Y] = meshgrid((1:size(obj.image,2))-floor(size(obj.image,2)/2)-1, (1:size(obj.image,1))-floor(size(obj.image,1)/2)-1);
                
                for ii = 1:C.NumObjects
                    
                    r(ii) = min(sqrt(X(C.PixelIdxList{ii}).^2 + Y(C.PixelIdxList{ii}).^2)); 
                    
                end
                
                [~,idx] = min(r); % find the connected component that is nearest to the center of the cutout
                
                N = numel(C.PixelIdxList{idx}); 
                
            end
            
            obj.num_pixels = N;
            
            obj.fwhm = util.img.fwhm(obj.image, 'method', 'filters'); 
            
            if obj.num_pixels==1
                obj.type = 'bad pixel';
            elseif obj.num_frames==1
                obj.type = 'cosmic ray';
            else
                obj.type = 'flare';
            end
            
        end
        
        function flux_norm = getNormFlux(obj)
            
            indices = obj.frame_index;
            for ii = 1:obj.interval_peak
                
                idx = obj.frame_index - ii;
                if idx<1, break; end
                
                indices = [indices idx]; 
                    
            end
            
            for ii = 1:obj.interval_peak
                
                idx = obj.frame_index + ii;
                if idx>length(obj.flux)
                    break;
                end
                
                indices = [indices idx]; 
                    
            end
            
            indices = sort(indices);
            flux_bg = obj.flux;
            flux_bg(indices) = NaN;
            
            obj.mean = nanmean(flux_bg);
            obj.std = nanstd(flux_bg); 
            
            flux_norm = (obj.flux-obj.mean)./obj.std;
            
        end
        
    end
    
    methods % other utilities (like save)
        
        function save(obj, filename)
            
            if nargin<2 || isempty(filename)
                error('Must supply a file name!');
            end
            
            event = obj;
            
            save(filename, 'event'); 
            
        end
        
        function saveDialog(obj, filename)
            
            if nargin<2 || isempty(filename)
                
                filename = 'flare'; 

                run_name = ''; 

                if ~isempty(obj.folder)
                    [a, b] = fileparts(obj.folder);
                    [~, c] = fileparts(a);
                    run_name = [c '_' b];
                end
                
                if ~isempty(run_name)
                    filename = [filename '_' run_name];
                end
                
                filename = sprintf('%s_id%03d', filename, obj.serial);
                
            end
            
            [filepath,filename,ext] = fileparts(filename);
            
            if ~isempty(ext)
                filename = [filename ext];
            else
                filename = [filename '.mat']; 
            end
            
            if isempty(filepath)
                filepath = fullfile(getenv('DATA'), 'WFAST/saved/flares/'); 
            end
            
            [filename, filepath] = uiputfile(fullfile(filepath, filename)); 
            
            if ~isequal(filename, 0)
                obj.save(fullfile(filepath, filename)); 
            end
            
        end
        
    end
    
    methods % plotting tools / GUI
        
        function show(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('index', []); 
            input.input_var('cutouts', 9, 'number'); 
            input.input_var('parent', []); 
            input.scan_vars(varargin{:}); 
            
            if isempty(input.parent)
                input.parent = gcf;
            end
            
            if ~isempty(input.index) % user input overrides all
                input.parent.UserData = input.index;
            elseif ~isempty(input.parent.UserData) % no user input, recover latest index
                input.index = input.parent.UserData;
                if input.index>length(obj)
                    input.index = 1;
                end
            else % the default is to start with the first index
                input.parent.UserData = 1;
                input.index = 1;
            end
            
            this = obj(input.index);
            
            clf(input.parent); 
            
            panel_filename = uipanel(input.parent, 'Units', 'Normalized', 'Position', [0 0.9 1 0.1]);
            
            uicontrol(panel_filename, 'Style', 'pushbutton', 'String', this.filename, ...
                'Units', 'Normalized', 'Position', [0 0 1 1], 'FontSize', 10); 
            
            panel_flux = uipanel(input.parent, 'Units', 'Normalized', 'Position', [0 0 0.5 0.9]); 
            
            uicontrol(panel_flux, 'Style', 'pushbutton', 'Units', 'Normalized', 'Position', [0.3 0.05 0.1 0.1], ...
                'String', '-', 'Callback', @obj.prevShow, 'FontSize', 16); 
            
            uicontrol(panel_flux, 'Style', 'edit', 'Units', 'Normalized', 'Position', [0.4 0.05 0.3 0.1], ...
                'String', sprintf('%d / %d ', input.index, length(obj)), 'Callback', @obj.chooseShow, 'FontSize', 16); 
            
            uicontrol(panel_flux, 'Style', 'pushbutton', 'Units', 'Normalized', 'Position', [0.7 0.05 0.1 0.1], ...
                'String', '+', 'Callback', @obj.nextShow, 'FontSize', 16); 
            
            uicontrol(panel_flux, 'Style', 'pushbutton', 'Units', 'Normalized', 'Position', [0.05 0.05 0.2 0.1], ...
                'String', 'save', 'Callback', @obj.saveCallback, 'FontSize', 16); 
            
            ax = axes('Parent', panel_flux, 'Position', [0.2 0.35 0.7 0.55]);
            
            plot(ax, this.flux, 'LineWidth', 2); 
            
            xlabel(ax, 'frame number'); 
            ylabel(ax, 'flux'); 
            title(ax, sprintf('peak = %4.2f | pos= %d %d | pix= %d | rms= %4.2f', this.peak, this.pos(1), this.pos(2), obj(input.index).num_pixels, sqrt(this.pixel_var))); 
            ax.FontSize = 14;
            
            panel_cutouts = uipanel(input.parent, 'Units', 'Normalized', 'Position', [0.5 0 0.5 0.9]); 
            
            ax = util.plot.show_cutouts(this.cutouts, 'parent', panel_cutouts, 'frame', this.frame_index, 'number', input.cutouts); 
            
            for ii = 1:length(ax)
                if ~isempty(ax{ii}.UserData) && ax{ii}.UserData==this.frame_index
                    
                    hold(ax{ii}, 'on'); 
                    x = size(this.image,2)/2 +0.5; % + this.offset(this.frame_index,1);
                    y = size(this.image,1)/2 +0.5; % + this.offset(this.frame_index,2);
                    plot(ax{ii}, x, y, 'ro', 'MarkerSize', 25); 
                    hold(ax{ii}, 'off'); 
                    
                end
            end
            
        end
        
        function nextShow(obj, hndl, ~)
            
            idx = hndl.Parent.Parent.UserData;
            
            idx = idx + 1;
            if idx>length(obj)
                idx = 1;
            end
            
            obj.show('parent', hndl.Parent.Parent, 'index', idx); 
            
        end
        
        function prevShow(obj, hndl, ~)
            
            idx = hndl.Parent.Parent.UserData;
            
            idx = idx - 1;
            if idx<1
                idx = length(obj);
            end
            
            obj.show('parent', hndl.Parent.Parent, 'index', idx); 
            
        end
        
        function chooseShow(obj, hndl, ~)
            
            idx = util.text.extract_numbers(hndl.String);
            
            if ~isempty(idx)
                idx = idx{1}; 
            end
            
            if ~isempty(idx)
                idx = idx(1);
            end
            
            if isnumeric(idx) && ~isempty(idx)
                obj.show('parent', hndl.Parent.Parent, 'index', idx); 
            end
            
        end
        
        function saveCallback(obj, hndl, ~)
            
            idx = hndl.Parent.Parent.UserData;
            
            obj(idx).saveDialog;
            
        end
        
    end    
    
end

