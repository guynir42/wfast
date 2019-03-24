classdef Lightcurves < handle

    properties(Transient=true)
        
    end
    
    properties % objects
        
    end
    
    properties % inputs/outputs
        
        frame_index;
        num_frames;
        
    end
    
    properties % switches/controls
        
        show_what = 'raw'; % can choose "raw" or "cal" or "both".
        use_show_all = 0; % if turned on, will show entire run, if off will show just what was already recorded. 
        show_num_stars = 10; % up to this number of stars are shown.
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
                
        fluxes;
        fluxes_cal;
        weights;
        offsets_x;
        offsets_y;
        widths;
        backgrounds;
        
    end
    
    properties(Hidden=true)
            
        fluxes_full;
        fluxes_cal_full;
        weights_full;
        offsets_x_full;
        offsets_y_full;
        widths_full;
        backgrounds_full;
        
        show_what_list = {'raw', 'cal', 'both'};
        
        version = 1.00;
        
    end
    
    methods % constructor
        
        function obj = Lightcurves(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'Lightcurves')
                if obj.debug_bit, fprintf('Lightcurves copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit, fprintf('Lightcurves constructor v%4.2f\n', obj.version); end
            
            end
            
        end
        
    end
    
    methods % reset/clear
        
        function reset(obj)
            
            obj.fluxes_full = [];
            obj.fluxes_cal_full = [];
            obj.weights_full = [];
            obj.offsets_x_full = [];
            obj.offsets_y_full = [];
            obj.widths_full = [];
            obj.backgrounds_full = [];
            
            obj.frame_index = 0;
            
        end
        
    end
    
    methods % getters
        
        function val = get.fluxes(obj)
            
            val = obj.fluxes_full(1:obj.frame_index,:);
            
        end
        
        function val = get.fluxes_cal(obj)
            
            val = obj.fluxes_cal_full(1:obj.frame_index,:);
            
        end
        
        function val = get.weights(obj)
            
            val = obj.weights_full(1:obj.frame_index,:);
            
        end
        
        function val = get.offsets_x(obj)
            
            val = obj.offsets_x_full(1:obj.frame_index,:);
            
        end
        
        function val = get.offsets_y(obj)
            
            val = obj.offsets_y_full(1:obj.frame_index,:);
            
        end
        
        function val = get.widths(obj)
            
            val = obj.widths_full(1:obj.frame_index,:);
            
        end
        
        function val = get.backgrounds(obj)
            
            val = obj.backgrounds_full(1:obj.frame_index,:);
            
        end
        
    end
    
    methods % setters
        
    end
    
    methods % calculations
        
        function startup(obj, num_frames, num_stars)
            
            obj.fluxes_full = NaN(num_frames, num_stars);
            obj.fluxes_cal_full = NaN(num_frames, num_stars);
            obj.weights_full = NaN(num_frames, num_stars);
            obj.offsets_x_full = NaN(num_frames, num_stars);
            obj.offsets_y_full = NaN(num_frames, num_stars);
            obj.widths_full = NaN(num_frames, num_stars);
            obj.backgrounds_full = NaN(num_frames, num_stars);
            
        end
        
        function input(obj, varargin)
            
            input = util.text.InputVars;
            input.use_ordered_numeric = 1;
            input.input_var('fluxes', [], 'fluxes_raw');
            input.input_var('weights', []);
            input.input_var('offsets_x', []);
            input.input_var('offsets_y', []);
            input.input_var('widths', []);
            input.input_var('backgrounds', []);
            input.scan_vars(varargin{:});
            
            N = size(input.fluxes,1);
            
            obj.fluxes_full(obj.frame_index+1:obj.frame_index+N,:) = input.fluxes;
            
            if ~isempty(input.weights)
                obj.weights_full(obj.frame_index+1:obj.frame_index+N,:) = input.weights;
            end
            
            if ~isempty(input.offsets_x)
                obj.offsets_x_full(obj.frame_index+1:obj.frame_index+N,:) = input.offsets_x;
            end
            
            if ~isempty(input.offsets_y)
                obj.offsets_y_full(obj.frame_index+1:obj.frame_index+N,:) = input.offsets_y;
            end
            
            if ~isempty(input.widths)
                obj.widths_full(obj.frame_index+1:obj.frame_index+N,:) = input.widths;
            end
            
            if ~isempty(input.backgrounds)
                obj.backgrounds_full(obj.frame_index+1:obj.frame_index+N,:) = input.backgrounds;
            end
            
            obj.frame_index = obj.frame_index + N;
            
        end
        
    end
    
    methods % plotting tools / GUI
        
    end    
    
end

 