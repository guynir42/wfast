classdef CutoutStorage < handle
% Keep a bunch of cutouts for different stars, across a full run, along
% with the results from the regular photometric pipeline. 
% This is useful for trying to compare photometric methods to the existing
% code. E.g., training a deep-network to get more stable photometry. 
% This object is given to the img.Analysis object which calls the input()
% method, giving it the img.Photomety object as an input. 
% To enable this, set the Analysis object's use_cutouts_store=1.

    properties(Transient=true)
        
    end
    
    properties % objects
        
        head@head.Header;
        this_input@util.text.InputVars; 
        
    end
    
    properties % inputs/outputs
        
        cutouts; % calibrated cutouts for this batch
        flux; % raw flux for this batch
        flux_rem_bg; % flux after removing background,
        auxiliary; % other photometric measurements
        juldates; % use julian date of flux measurements to get airmass
        
        cutouts_all; % storage across the entire run
        flux_all; % flux for entire run, after any corrections e.g., zero-point
        flux_mean; % mean flux for each batch
        flux_corr; % after applying all sorts of corrections
        aux_all; % auxiliary measurements for entire run
        jul_all; % julian date for all flux measurements
        
        frame_counter = 0;
        
    end
    
    properties % switches/controls
        
        aux_names = {'errors', 'areas', 'backgrounds', 'variances', 'offsets_x', 'offsets_y', 'centroids_x', 'centroids_y', 'widths', 'bad_pixels', 'flags'}; % add more aux measurements if you want! 
        aux_indices; % struct with field=number for each of the above aux names
        
        frame_indices = [1 50]; % choose which frames in each batch get saved
        star_indices = []; % empty vector means save all stars
        aperture_idx = 1; % which flux to take in case of multiple fluxes/apertures
        
        use_zero_point = 1; % also do some processing on the flux, using each frame's zero point 
        use_airmass = 1; % also make an adjustment based on airmass
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        
        
    end
    
    properties(Hidden=true)
       
        version = 1.00;
        
    end
    
    methods % constructor
        
        function obj = CutoutStorage(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'learn.CutoutStorage')
                if obj.debug_bit>1, fprintf('CutoutStorage copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit>1, fprintf('CutoutStorage constructor v%4.2f\n', obj.version); end
            
                obj.reset;
                
            end
            
        end
        
    end
    
    methods % reset/clear
        
        function reset(obj)
            
            obj.cutouts_all = [];
            obj.flux_all = [];
            obj.flux_mean = [];
            obj.flux_corr = []; 
            obj.aux_all = []; 
            obj.jul_all = [];
            
            obj.frame_counter = 0;
            
            obj.calcAuxIndices;
            
            obj.clear;
            
        end
        
        function clear(obj)
            
            obj.cutouts = [];
            obj.flux = [];
            obj.flux_rem_bg = []; 
            obj.auxiliary = [];
            obj.juldates = []; 
            
        end
        
        function calcAuxIndices(obj)
            
            obj.aux_indices = struct;
            
            for ii = 1:length(obj.aux_names)
                obj.aux_indices.(obj.aux_names{ii}) = ii; 
            end
            
        end
        
    end
    
    methods % getters
        
    end
    
    methods % setters
        
    end
    
    methods % calculations
        
        function val = makeInputVars(obj) % make an input object that can parse all photometric products and meta-data
            
            val = util.text.InputVars;
            val.use_ordered_numeric = 1; % can just give the numeric values in the correct order, without naming each one
            
            % these inputs are used for explicitely giving the photometric products
            val.input_var('fluxes', []);
            val.input_var('cutouts', []); 
            
            for ii = 1:length(obj.aux_names)
                if ~strcmp(obj.aux_names{ii}, 'fluxes')
                    val.input_var(obj.aux_names{ii}); 
                end
            end
            
            val.input_var('juldates', [], 'julian_dates'); 
            val.input_var('pars_struct', []); 
            
        end
        
        function input(obj, varargin)
            
            if isempty(varargin)
                error('Must provide photometric products like timestamps, fluxes and so on'); 
            elseif isa(varargin{1}, 'util.text.InputVars')
                obj.this_input = util.oop.full_copy(varargin{1}); % just take the input object some one else provided
                % what happens to any other varargin pairs? 
            elseif isa(varargin{1}, 'img.Photometry') % giving a photometry object is the prefered method! 
                
                obj.this_input = obj.makeInputVars;
                
                list = obj.this_input.list_added_properties;

                for ii = 1:length(list)
                    obj.this_input.(list{ii}) = varargin{1}.(list{ii}); 
                end
               
            else % just make an input object and parse the varargin pairs
                obj.this_input = obj.makeInputVars;
                obj.this_input.scan_vars(varargin{:}); 
            end
            
            % make sure we got fluxes, timestamps and all the aux data
            if isempty(obj.this_input.fluxes) || ...
                    isempty(obj.this_input.errors) || isempty(obj.this_input.areas) ||...
                    isempty(obj.this_input.backgrounds) || isempty(obj.this_input.variances) || ...
                    isempty(obj.this_input.offsets_x) || isempty(obj.this_input.offsets_y) || ...
                    isempty(obj.this_input.centroids_x) || isempty(obj.this_input.centroids_y) || ...
                    isempty(obj.this_input.widths) || isempty(obj.this_input.bad_pixels) || ...
                    isempty(obj.this_input.flags)

                error('Some photometric products are empty!'); 

            end
            
            if ~isempty(obj.head) && ~isempty(obj.this_input.pars_struct)
                obj.head.PHOT_PARS = obj.this_input.pars_struct; % update the header with the photometric parameters used in the analysis
            end
            
            if ndims(obj.this_input.fluxes)>3
                error('This class cannot handle more than 3D fluxes!'); 
            end
            
            if ndims(obj.this_input.widths)>3
                error('This class cannot handle more than 3D auxiliary measurements!'); 
            end
            
            obj.clear; % get rid of the data from last batch
            
            obj.frame_counter = obj.frame_counter + size(obj.this_input.fluxes,1); % how many frames we processed
            
            obj.cutouts = obj.this_input.cutouts; 
            
            obj.flux = obj.this_input.fluxes(:,:,obj.aperture_idx); 
            obj.juldates = obj.this_input.juldates; 
            
            % store the new aux data
            obj.auxiliary = NaN(size(obj.this_input.errors,1), size(obj.this_input.widths,2), length(obj.aux_names), 'like', obj.this_input.widths); % preallocate
            
            for ii = 1:length(obj.aux_names)
                new_matrix = obj.this_input.(obj.aux_names{ii});
                obj.auxiliary(:,:,ii) = new_matrix(:,:,obj.aperture_idx); % allow multiple apertures (3D aux matrices from photometry) 
            end
            
            % now process the fluxes as best you can
            obj.flux_rem_bg = obj.flux - obj.auxiliary(:,:,obj.aux_indices.areas).*obj.auxiliary(:,:,obj.aux_indices.backgrounds); % subtract background
            
            f = obj.flux_rem_bg; 

            for ii = 1:2
                
                s = nanstd(f); 
                m = nanmean(f); 
                bad_idx = abs(f-m)./s > 5; % remove 5 sigma outliers
                
                f(bad_idx) = NaN; % replace outliers with NaNs that are not accounted in the averages
                
            end
            
            F = nanmean(f); % batch average flux
            
            % save the data into full-run storage
            idx_frame = obj.frame_indices;
            if isempty(idx_frame)
                idx_frame = 1:size(f,1); % empty indices means grab everything!
            end
            
            idx_star = obj.star_indices;
            if isempty(idx_star)
               idx_star = 1:size(f,2); % empty indices means grab everything!
            end
            
            obj.cutouts_all = cat(3, obj.cutouts_all, obj.cutouts(:,:,idx_frame,idx_star));
            obj.flux_all = cat(1, obj.flux_all, obj.flux_rem_bg(idx_frame, idx_star));
            obj.flux_mean = cat(1, obj.flux_mean, repmat(F(1,idx_star), [length(idx_frame), 1])); % repmat a mean flux for each frame index / cutout saved
            
            obj.aux_all = cat(1, obj.aux_all, obj.auxiliary(idx_frame, idx_star,:)); 
            obj.jul_all = cat(1, obj.jul_all, obj.juldates(idx_frame));
            
        end
        
        function [Xvalues, Yvalues, Xvalidation, Yvalidation] = getData(obj, num_validation) % output the data as a mapping between cutouts->photometric measurements
            
            if nargin<2 ||isempty(num_validation)
                num_validation = ceil(size(obj.flux_mean,1)*0.2); % 20% of the data is used as validation
            end
            
            S = size(obj.cutouts_all); 
            N = S(3).*S(4); 
            
            A = obj.aux_all(:,:,[obj.aux_indices.backgrounds, obj.aux_indices.offsets_x, obj.aux_indices.offsets_y, obj.aux_indices.widths]); 
            A = reshape(A, [N,size(A,3),1]);
            
%             F = oobj.flux_mean; 
            F = obj.flux_all./obj.flux_mean; 
            F = reshape(F, [N, 1]);
            
            Xvalues = reshape(obj.cutouts_all, [S(1), S(2), 1, N]); 
            Yvalues = horzcat(F, A); 
            
            bad_idx = any(isnan(Yvalues),2); 
            
            Xvalues(:,:,:,bad_idx) = [];
            Yvalues(bad_idx,:) = []; 
            
            idx = randperm(size(Yvalues,1), num_validation); 
            
            Xvalidation = Xvalues(:,:,:,idx);
            Yvalidation = Yvalues(idx,:); 
            
            Xvalues(:,:,:,idx) = [];
            Yvalues(idx,:) = []; 
            
        end
        
    end
    
    methods % plotting tools / GUI
        
    end    
    
end

