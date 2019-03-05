classdef Manager < handle

    properties(Transient=true)
        
        gui;
        
    end
    
    properties % objects
        
        checker@obs.StatusChecker;
        
        dome;
        mount;
        
        weather;
        wind;
        humidity;
        temperature;
        
    end
    
    properties % inputs/outputs
        
    end
    
    properties % switches/controls
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        
        
    end
    
    properties(Hidden=true)
       
        version = 1.00;
        
    end
    
    methods % constructor
        
        function obj = Manager(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'obs.Manager')
                if obj.debug_bit, fprintf('Manager copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit, fprintf('Manager constructor v%4.2f\n', obj.version); end
            
                obj.connect;
                obj.checker = obs.StatusChecker(obj);
                
            end
            
        end
        
        function connect(obj)
            
            try 
                obj.dome = obs.dome.AstroHaven;
            catch ME
                
                warning(ME.getReport);
                
                disp('Cannot connect to AstroHaven dome. Using simulator instead...');
                
                try
                    obj.dome = obs.dome.Simulator;
                catch ME
                    warning(ME.getReport);
                end
                
            end
            
            try 
                obj.mount = obs.mount.ASA;
            catch ME
                
                warning(ME.getReport);
                
                disp('Cannot connect to ASA mount. Using simulator instead...');
                
                try 
                    obj.mount = obs.mount.Simulator;
                catch ME
                    warning(ME.getReport);
                end
                
            end
            
            try 
                obj.weather = sens.Boltwood;
            catch ME
                
                warning(ME.getReport);
                
                disp('Cannot connect to Boltwood weather station. Using simulator instead...');
                
                try
                    obj.weather = obs.sens.Simulator;
                catch ME
                    warning(ME.getReport);
                end
                
            end
            
        end
        
    end
    
    methods % reset/clear
        
    end
    
    methods % getters
        
    end
    
    methods % setters
        
    end
    
    methods % calculations
        
    end
    
    methods % plotting tools / GUI
        
    end    
    
end

