classdef Manager < handle

    properties(Transient=true)
        
        gui;
        
    end
    
    properties % objects
        
        log@util.sys.Logger;
        
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
                
                obj.log = util.sys.Logger('Top_level_manager'); % keep track of commands given and errors received... 
                
                obj.connect; % connect to all hardware
                
                obj.checker = obs.StatusChecker(obj); % has timers that check weather/hardware status
                
            end
            
        end
        
        function connect(obj)
            
            obj.connectDome;
%             obj.connectMount;
            obj.connectBoltwood;
            
        end
        
        function connectDome(obj)
            
            obj.log.input('Connecting to dome.');
            
            try 
                obj.dome = obs.dome.AstroHaven;
            catch ME
                
                obj.log.error(ME.getReport);
                obj.log.input('Connecting dome simulator.');
                
                warning(ME.getReport);
                
                disp('Cannot connect to AstroHaven dome. Using simulator instead...');
                
                try
                    obj.dome = obs.dome.Simulator;
                catch ME
                    obj.log.error(ME.getReport);
                    warning(ME.getReport);
                end
                
            end
            
        end
        
        function connectMount(obj)
            
            obj.log.input('Connecting to mount.');
            
            try 
                obj.mount = obs.mount.ASA;
            catch ME
                
                obj.log.error(ME.getReport);
                obj.log.input('Connecting to mount simulator.');
                
                warning(ME.getReport);
                
                disp('Cannot connect to ASA mount. Using simulator instead...');
                
                try 
                    obj.mount = obs.mount.Simulator;
                catch ME
                    obj.log.error(ME.getReport);
                    warning(ME.getReport);
                end
                
            end
            
        end
        
        function connectBoltwood(obj)
            
            obj.log.input('Connecting to Boltwood weather station.');
            
            try 
                obj.weather = obs.sens.Boltwood;
            catch ME
                
                obj.log.error(ME.getReport);
                obj.log.input('Connecting to weather simulator.');
                warning(ME.getReport);
                
                disp('Cannot connect to Boltwood weather station. Using simulator instead...');
                
                try
                    obj.weather = obs.sens.Simulator;
                catch ME
                    obj.log.error(ME.getReport);
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

