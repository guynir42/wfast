classdef StarHours < handle

    properties(Transient=true)
        
    end
    
    properties % objects
        
    end
    
    properties % inputs/outputs
        
        snr_bin_edges; % this is calculated from the parameters below
        star_bin_edges; % this is calculated after we know how many stars we've got
        
        runtime = 0; % total runtime since last reset()
        histogram; % the number of useful seconds accumulated for each star and each S/N bin (this is saved after subtracting losses)
        losses_exclusive; % same as histogram, only counting the losses due to each cut (exclusive means times where ONLY this cut was responsible)
        losses_inclusive; % same as histogram, only counting the losses due to each cut (inclusive means times where this cut was ALSO responsible)
        
        cut_names = {}; % get these from the QualityChecker
        cut_indices = {}; % get these from the QualityChecker
        
        % these are temporary values used each batch:
        timestamps; % timestamps for the extended region
        flux; % flux for the extended region
        idx_start; % start of the search region
        idx_end; % end of the search region
        star_snr; % get the S/N for each star based on the extended region
        time_step; % average dt (0.04s)
        batch_length; % how many seconds are in the batch, total
        flags; % a matrix with all the regions flagged for each cut
        
    end
    
    properties % switches/controls
        
        snr_bin_min = 5; % stars with lower S/N are not even tested for events (here S/N is calculated on the raw fluxes)
        snr_bin_width = 0.5; % for the star-hour histogram
        snr_bin_max = 100; % biggest value we expect to see in the star-hour histogram
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        
        
    end
    
    properties(Hidden=true)
       
        star_snr_dist; % output of histcounts2(). Should be 1 in the right star and S/N bin, and zeros elsewhere. 
        
        version = 1.00;
        
    end
    
    methods % constructor
        
        function obj = StarHours(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'trig.StarHours')
                if obj.debug_bit>1, fprintf('StarHours copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit>1, fprintf('StarHours constructor v%4.2f\n', obj.version); end
            
                obj.reset;
                
            end
            
        end
        
    end
    
    methods % reset/clear
        
        function reset(obj)
            
            obj.runtime = 0;
            obj.histogram = [];
            obj.losses_exclusive = [];
            obj.losses_inclusive = [];
            obj.snr_bin_edges = [];
            obj.star_bin_edges = []; 
            
            obj.clear; 
            
        end
        
        function clear(obj)
            
            obj.timestamps = [];
            obj.flux = [];
            obj.idx_start = [];
            obj.idx_end = [];
            obj.star_snr = [];
            obj.time_step = []; 
            obj.batch_length = [];
            obj.flags = []; 
            obj.star_snr_dist = [];
            
        end
        
    end
    
    methods % getters
        
    end
    
    methods % setters
        
    end
    
    methods % calculations
        
        function input(obj, checker)
            
            obj.clear;
            
            % get the data from checker
            obj.timestamps = checker.extended_timestamps;
            obj.flux = checker.extended_flux;
            obj.idx_start = checker.search_start_idx; 
            obj.idx_end = checker.search_end_idx; 
        
            obj.flags = checker.cut_flag_matrix; 
            
            obj.cut_names = checker.cut_names; 
            obj.cut_indices = checker.cut_indices;
        
            % calculate some stuff based on that
            obj.star_snr = nanmean(obj.flux,1)./nanstd(obj.flux); % S/N on the extended region
            obj.time_step = median(diff(obj.timestamps)); 
            obj.batch_length = (obj.idx_end - obj.idx_start + 1).*obj.time_step; 
            
            if isempty(obj.star_bin_edges) || isempty(obj.snr_bin_edges) || isempty(obj.histogram)
                obj.makeHistograms;
            end
            
            % fill up the cumulative stuff for the whole run:
            obj.runtime = obj.runtime + obj.batch_length; 
            
            obj.star_snr_dist = histcounts2(obj.star_snr, obj.star_bin_edges(1:end-1), obj.snr_bin_edges, obj.star_bin_edges-0.5); % resulting histogram is stars on the X axis and S/N on the Y axis
            bad_times = checker.bad_times; % this has "true" in every frame/star that is disqualified for anything
            lost_time = sum(bad_times).*obj.time_step; 
            
            obj.histogram = obj.histogram + obj.star_snr_dist.*(obj.batch_length - lost_time); % add the amount of time for each star in the right bin
            
%             flag_also = false(size(bad_times,1), size(bad_times,2), length(obj.cut_names)); % logical "true" in each frame/star where this cut is also or only included
%             flag_only = false(size(bad_times,1), size(bad_times,2), length(obj.cut_names)); % logical "true" in each frame/star where this cut is only included (with no others!)
            
            flag_only = false(size(obj.flags)); % logical "true" in each frame/star where this cut is only included (with no others!)

            for ii = 1:length(obj.cut_names)
                
                num_xors = sum(xor(obj.flags(:,:,ii), obj.flags),3); % each frame/star counts the number of places in other cuts that were not "true" at the same time as the current cut
                flag_only(:,:,ii) = num_xors==length(obj.cut_names)-1; % for each frame/star if all other cuts have zero here, then the xor would be "true" for all other cuts except the current cut
            
                obj.losses_inclusive(:,:,ii) = obj.losses_inclusive(:,:,ii) + obj.star_snr_dist.*sum(obj.flags(:,:,ii)).*obj.time_step; % add the amount of time (per star, per S/N) for this cut (even if this time is shared with other cuts)
                obj.losses_exclusive(:,:,ii) = obj.losses_exclusive(:,:,ii) + obj.star_snr_dist.*sum(flag_only(:,:,ii)).*obj.time_step; % add the amount of time (per star, per S/N) for this cut (even if this time is shared with other cuts)
                
            end
            
        end
        
        function makeHistograms(obj)
            
            obj.star_bin_edges = 1:size(obj.flux,2)+1;
            obj.snr_bin_edges = obj.snr_bin_min:obj.snr_bin_width:obj.snr_bin_max;
            
            obj.histogram = zeros(length(obj.snr_bin_edges)-1, length(obj.star_bin_edges)-1); 
            obj.losses_exclusive = zeros(length(obj.snr_bin_edges)-1, length(obj.star_bin_edges)-1, length(obj.cut_names));
            obj.losses_inclusive = zeros(length(obj.snr_bin_edges)-1, length(obj.star_bin_edges)-1, length(obj.cut_names)); 
            
        end
        
    end
    
    methods % plotting tools / GUI
        
        function showHistogram(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('type', []); % choose cut type or leave empty to get the good hours
            input.input_var('exclusive', false); % if false, show inclusive losses
            input.input_var('stars', []); 
            input.input_var('max_snr', obj.snr_bin_max); 
            input.input_var('min_snr', obj.snr_bin_min); 
            input.input_var('rebin', []); 
            input.input_var('sum', false);
            input.input_var('hours', true); 
            input.input_var('axes', [], 'axis'); 
            input.input_var('font_size', 18); 
            input.input_var('log', true, 'use_log'); 
            input.input_var('color', []); 
            input.input_var('alpha', 0.7); 
            input.scan_vars(varargin{:}); 
            
            if isempty(input.axes)
                input.axes = gca;
            end
            
            if isempty(input.type)
                
                values = obj.histogram;
                str = 'Effective star-time'; 
                
            else
                if ischar(input.type)
                    input.type = obj.cut_indices.(input.type); 
                end

                if input.exclusive
                    values = obj.losses_exclusive(:,:,input.type); 
                    exc = '(exclusively)';
                    
                else
                    values = obj.losses_inclusive(:,:,input.type); 
                    exc = '(inclusively)'; 
                end
                
                str = sprintf('Time lost %s to "%s"', exc, strrep(obj.cut_names{input.type}, '_', ' '));
                
            end
            
            units = '[seconds]'; 
            
            if input.hours
                values = values./3600; % transform to hours instead of seconds!
                units = '[hours]'; 
            end
            
            if input.sum
                
                values = sum(values,2);
                h = bar(input.axes, obj.snr_bin_edges(1:end-1), values, 'FaceAlpha',input.alpha); 
                
                if ~isempty(input.color)
                    h.FaceColor = input.color;
                end
                
                mx = max(values); 
                
                if mx<=0
                    mx = 1;
                else
                    mx = mx.*1.1;
                end
                
                input.axes.YLim = [0 mx]; 
                
                if input.log
                    input.axes.YScale = 'log'; 
                end
                
                h.DisplayName = [str ' ' units];
            
                legend(input.axes, 'Location', 'NorthEast'); 
            
            else
                
                cla(input.axes); 
                
                imagesc(input.axes, 'XData', obj.snr_bin_edges(1:end-1)+obj.snr_bin_width/2, 'CData', values');
                ylabel(input.axes, 'Star index');                 
                colorbar(input.axes); 
                input.axes.YLim = [0 length(obj.star_bin_edges)-1]; 
                
                input.axes.YScale = 'linear'; 
                
                if input.log
                    input.axes.ColorScale = 'log'; 
                end
                
                legend(input.axes, 'off'); 
              
                util.plot.inner_title([str ' ' units], 'ax', input.axes, 'Position', 'NorthEast', 'Color', 'white'); 
            
            end
            
            input.axes.XLim = [input.min_snr, input.max_snr-obj.snr_bin_width];
            xlabel(input.axes, 'S/N'); 
            
            input.axes.FontSize = input.font_size;
            axis(input.axes, 'normal'); 
            axis(input.axes, 'xy'); 
            
        end
        
    end
    
end







