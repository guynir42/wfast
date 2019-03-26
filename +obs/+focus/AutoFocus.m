classdef AutoFocus < handle

    properties(Transient=true)
        
        fig; % for displaying the focus curve
        ax; 
        
    end
    
    properties % objects
        
    end
    
    properties % inputs/outputs
        
        pos; % focuser position at each iteration
        widths; % one or multiple measurements per focus position
        weights; % weight of each sampling point (proportional to star flux)
        xy_pos; % position of sampling points
        
        fit_results = {};
        min_positions; % minimum of curves for each star
        
        surface_coeffs;
        
        found_pos; % put your focuser to this position
        found_tip; % put your focuser to this tip value
        found_tilt; % put your focuser to this tilt value
        
    end
    
    properties % switches/controls
        
        step = 0.01;
        range = 0.1;
        
        angle = 0; % between tip axis and pixel y axis (degrees)
        spider_diameter = 100; % in cm
        pixel_size = 6.5; % in microns
        num_pixels = 2000; % across the sensor (roughly)
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        
        
    end
    
    properties(Hidden=true)
       
        version = 1.00;
        
    end
    
    methods % constructor
        
        function obj = AutoFocus(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'obs.focus.AutoFocus')
                if obj.debug_bit, fprintf('AutoFocus copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit, fprintf('AutoFocus constructor v%4.2f\n', obj.version); end
            
            end
            
        end
        
    end
    
    methods % reset/clear
        
        function reset(obj)
            
            obj.pos = [];
            obj.widths = [];
            obj.weights = [];
            obj.xy_pos = [];
            
            obj.fit_results = {};
            obj.min_positions = [];
            obj.surface_coeffs = [];
            
            obj.found_pos = [];
            obj.found_tip = [];
            obj.found_tilt = [];
            
            obj.clear;
            
        end
        
        function clear(obj)
            
        end
        
    end
    
    methods % getters
        
    end
    
    methods % setters
        
    end
    
    methods % calculations
        
        function calculate(obj)
            
            obj.fitCurves;
%             obj.fitSurface;
%             obj.findPosTipTilt;
            obj.findPosOnly;
            
        end
        
        function fitCurves(obj)
            
            obj.fit_results = {};
            
            for ii = 1:size(obj.widths, 2)
                
                obj.fit_results{ii} = fit(obj.pos, obj.widths(:,ii), 'poly2');
                
                obj.min_positions(ii) = -obj.fit_results{ii}.p2./(2*obj.fit_results{ii}.p1);
                
            end
            
        end
        
        function fitSurface(obj)
            
            m = size(obj.xy_pos,1); % number of measurements
                        
            B = obj.min_positions; % measured best position for each location
            w = obj.weights; % flux of each star gives the level of conifidence
            
            xy_rot = obj.xy_pos*[cosd(obj.angle), -sind(obj.angle); sind(obj.angle), cosd(obj.angle)]; 
                        
            A = [ones(m,1) xy_rot(:,1) xy_rot(:,2)]; % design matrix! 

            obj.surface_coeffs = lscov(A,B,w); % coeffs are: piston, x (tilt) and y (tip)
            
        end
        
        function findPosOnly(obj)
            
            obj.found_pos = mean(obj.min_positions);
        
        end
        
        function findPosTipTilt(obj)
            
            obj.found_pos = obj.surface_coeffs(1);
            obj.found_pos = mean(obj.min_positions);
            if obj.debug_bit, fprintf('BEST POS: mean= %f | surface piston term= %f\n', obj.found_pos, mean(obj.min_positions)); end
            
            obj.found_tilt = obj.surface_coeffs(2).*obj.spider_diameter.*1e4./obj.pixel_size;
            obj.found_tip = obj.surface_coeffs(3).*obj.spider_diameter.*1e4./obj.pixel_size;
            
            if obj.debug_bit, fprintf('BEST TIP= %f | BEST tilt= %f\n', obj.found_tip, obj.found_tilt); end
            
        end
            
    end
    
    methods % plotting tools / GUI
        
        function plot(obj, varargin)
            
            if isempty(obj.fig) || ~isvalid(obj.fig)
                obj.fig = figure;
            end
            
            if isempty(obj.ax) || ~isvalid(obj.ax)
                obj.ax = axes('Parent', obj.fig);
            end
            
            plot(obj.ax, obj.pos, obj.widths);
            
            hold(obj.ax, 'on');
            
            for ii = 1:length(obj.fit_results)
                
                plot(obj.ax, obj.pos, feval(obj.fit_results{ii}, obj.pos));
                
            end
            
            hold(obj.ax, 'off');
            
            
        end
        
    end
    
end

