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
        
        raw;
        cal;
        
    end
    
    properties(Hidden=true)
       
        raw_full;
        cal_full;
        has_cal = 0;
        
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
            
            obj.raw_full = [];
            obj.cal_full = [];
            obj.has_cal = 0;
            obj.frame_index = 0;
            
        end
        
    end
    
    methods % getters
        
        function val = get.raw(obj)
            
            val = obj.raw_full(1:obj.frame_index,:);
            
        end
        
        function val = get.cal(obj)
            
            val = obj.cal_full(1:obj.frame_index,:);
            
        end
        
    end
    
    methods % setters
        
    end
    
    methods % calculations
        
        function startup(obj, num_frames, num_stars)
            
            obj.raw_full = NaN(num_frames, num_stars);
            obj.cal_full = NaN(num_frames, num_stars);
            
        end
        
        function input(obj, fluxes_raw, fluxes_cal)
            
            if nargin<3 || isempty(fluxes_cal)
                fluxes_cal = [];
            end
            
            N = size(fluxes_raw,1);
            
            obj.raw_full(obj.frame_index+1:obj.frame_index+N,:) = fluxes_raw;
            
            
            if ~isempty(fluxes_cal)
                
                obj.cal_full(obj.frame_index+1:obj.frame_index+N,:) = fluxes_cal;
                obj.has_cal = 1;
                
            end
                
            obj.frame_index = obj.frame_index + N;
            
        end
        
    end
    
    methods % plotting tools / GUI
        
    end    
    
end

 