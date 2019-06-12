classdef ASA < handle

    properties(Transient=true)
        
    end
    
    properties % objects
        
        hndl;
        
        target@head.Ephemeris;
        
        ard@obs.sens.ScopeAssistant;
        
        log@util.sys.Logger;
        
    end
    
    properties % inputs/outputs
        
        status = 0;
        
    end
    
    properties % switches/controls
        
        limit_alt = 15; % degrees
        
        use_accelerometer = 1;
        use_ultrasonic = 0;
        
        step_arcsec = 5;
        
        brake_bit = 1; 
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        RA_target; % in degrees        
        DE_target; % in degrees
        
        RA_target_hex;
        DE_target_hex;
        
        HA_target;
        ALT_target;
        AZ_target;
        
        RA;
        DE;        
        
        RA_hex;
        DE_hex;
        
        HA;
        ALT;
        AZ;
        
        LST;
        LST_hex;
        
        % can we get motor status from hndl??
        tracking;
        
        rate;
        
    end
    
    properties(Hidden=true)
       
        version = 1.01;
        
    end
    
    methods % constructor
        
        function obj = ASA(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'obs.mount.ASA')
                if obj.debug_bit, fprintf('ASA copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit, fprintf('ASA constructor v%4.2f\n', obj.version); end
            
            end
            
            obj.log = util.sys.Logger('ASA_mount', obj);
            
            try
            
                obj.target = head.Ephemeris;

                obj.connect;

                obj.update;

            catch ME
                obj.log.error(ME.getReport);
                rethrow(ME);
            end
            
        end
        
        function delete(obj)
           
            obj.disconnect;
            
        end
        
        function connect(obj)
            
            obj.log.input('Connecting to mount');
            
            try 
                
                if ~isempty(obj.hndl)
                    obj.disconnect;
                    pause(0.1);
                end
                
                obj.loadServer;
                
                obj.hndl = actxserver('AstrooptikServer.Telescope');
            
                obj.hndl.SiteLatitude  = obj.target.latitude;
                obj.hndl.SiteLongitude = obj.target.longitude;
                
                obj.hndl.Connected = 1;
            
                obj.hndl.MotorOn;
                
                if obj.use_accelerometer
                    
                    if isempty(obj.ard) 
                        try
                            obj.ard = obs.sens.ScopeAssistant;
                            obj.ard.telescope = obj;
                        catch ME
                            obj.use_accelerometer = 0;
                            warning(ME.getReport);
                        end
                    end
                end
                
                if obj.use_accelerometer
                    obj.ard.connect;
                end
                
            catch ME
                obj.log.error(ME.getReport);
                rethrow(ME);
            end
            
        end
        
        function loadServer(obj)
            
%             system('C:\Program Files (x86)\Autoslew\AstroOptikServer.exe &'); % call the system command to load the server outside of matlab 
            system('D:\matlab\wfast\+obs\+mount\launch_server.bat &'); % call the batch file to load the server then exit the cmd
            
            tic;
            
            for ii = 1:100
                
                [~, str] = system('tasklist /fi "imagename eq AstroOptikServer.exe"'); % check 
                
                str = strip(str);
                idx = strfind(str, 'AstroOptikServer.exe');
                
                if isempty(idx)
                    continue; % can't find the server on the tasklist
                else
                    num = util.text.extract_numbers(str(idx:end)); % found the server, now check the memory usage
                    if num{1}(3)>30
                        return; % server has loaded and is connected to telescope
                    end
                    % if this last condition fails, repeatedly, it means
                    % the server is stuck on a dialog... 
                end
                
                pause(0.05);
            
            end
            
            obj.killServer;
            error('Timeout while waiting for AstroOptikServer.exe to load for %f seconds!', toc); 
            
        end
        
        function disconnect(obj)
            
            obj.log.input('Disconnecting from mount');
            
            try
                obj.killServer;
                delete(obj.hndl);
                obj.hndl = [];
            catch ME
                obj.log.error(ME.getReport);
                rethrow(ME);
            end
            
            
        end
        
        function killServer(obj)
            
            [ret, str] = system('taskkill /fi "imagename eq AstroOptikServer.exe" /f'); % force killing of the server.
            
        end
        
    end
    
    methods % reset/clear
        
        function reset(obj)
            
            obj.target.reset;
            
        end
        
    end
    
    methods % getters
        
        function val = get.RA_target(obj)
            % get RA of user target in deg
            val = rad2deg(obj.target.RA_rad);
            
        end
        
        function val = get.DE_target(obj)
            % get Dec of user target in deg
            val = rad2deg(obj.target.DE_rad);
            
        end
        
        function val = get.RA_target_hex(obj)
            
            val = obj.target.RA; 
            
        end
        
        function val = get.DE_target_hex(obj)
            
            val = obj.target.DE;
            
        end
        
        function val = get.HA_target(obj)
            
            val = obj.target.HA;
            
        end
        
        function val = get.ALT_target(obj)
            
            val = rad2deg(obj.target.ALT);
            
        end
        
        function val = get.AZ_target(obj)
            
            val = rad2deg(obj.target.AZ);
            
        end
        
        function val = get.RA(obj)
            
            val = obj.hndl.RightAscension/24*360;
            
        end
        
        function val = get.DE(obj)
            
            val = obj.hndl.Declination;
            
        end
        
        function val = get.RA_hex(obj)
            
            val = head.Ephemeris.rad2ra(deg2rad(obj.RA));
            
        end
        
        function val = get.DE_hex(obj)
            
            val = head.Ephemeris.rad2dec(deg2rad(obj.DE));
            
        end
        
        function val = get.ALT(obj)
            
            val = obj.hndl.Altitude;
            
        end
        
        function val = get.AZ(obj)
            
            val = obj.hndl.Azimuth;
            
        end
        
        function val = get.HA(obj)
            
            val = obj.LST - obj.RA;
            
        end
        
        function val = get.LST(obj)
            
            val = obj.hndl.SiderealTime/24*360; % in degrees...
            
        end
        
        function val = get.LST_hex(obj)
            
            val = head.Ephemeris.rad2ra(deg2rad(obj.LST));
            
        end
        
        function val = latitutde(obj)
            
%             val = obj.hndl.SiteLatitude;
            val = obj.target.latitude;
            
        end
        
        function val = longitude(obj)
           
%             val = obj.hndl.SiteLongitude;
            val = obj.target.longitude;

        end
        
    end
    
    methods % setters
        
        function set.RA_target(obj, val)
            
            obj.hndl.TargetRightAscension = val*24/360;
            obj.target.RA = deg2rad(val);
%             obj.RA_target = val; % should we make this dependent?
            
        end
        
        function set.DE_target(obj, val)
            
            obj.hndl.TargetDeclination = val;
            obj.target.DE = deg2rad(val);
%             obj.DE_target = val; % should we make this dependent?
            
        end
        
    end
    
    methods % calculations / commands
        
        function inputTarget(obj, varargin)
            
            obj.target.input(varargin{:}); % Ephemeris now uses Eran's name resolver
            
        end
        
        function val = check_before_slew(obj)
            
            val = 0;
            
            obj.target.update;
            
            if obj.target.Alt_deg<obj.limit_alt
                return;
            end
            
            val = 1;
            
        end
        
        function val = check_while_moving(obj)
            
            val = 0;
            
            if obj.ALT<obj.limit_alt
                return;
            end
            
            val = 1;
            
        end
        
        function slew(obj, varargin)
            
            if isempty(obj.target.RA) || isempty(obj.target.DE) 
                error('Please provide a target with viable RA/DE');
            end
            
            obj.log.input(sprintf('Slewing to target. RA= %s | DE= %s | ALT= %4.2f', obj.target.RA, obj.target.DE, obj.target.ALT_deg));
            
            try 
                
                if ~obj.check_before_slew
                    error('Prechecks failed, aborting slew');
                end
                
                ra_hours_Jnow = obj.RA_deg_now./15; % convert to hours! 
                dec_deg_Jnow = obj.target.Dec_deg_now;
                
                obj.brake_bit = 0;
                obj.hndl.SlewToCoordinatesAsync(ra_hours_Jnow, dec_deg_Jnow);
                
                for ii = 1:100000
                    
                    pause(0.01);
                    
                    if obj.brake_bit
                        obj.stop;
                        break;
                    end
                    
                    if ~obj.check_while_moving
                        obj.emergency_stop;
                        break;
                    end
                    
                    if obj.hndl.Slewing==0
                        break;
                    end
                    
                end
            
                obj.tracking = 1;
                
            catch ME
                obj.log.error(ME.getReport);
                rethrow(ME);
            end
            
        end
        
        function adjustPosition(obj, RA_deg, DE_deg)
            
            error('This doesnt work, dont use it');
            
            obj.log.input(sprintf('Adjusting position by RA: %f deg | DE: %f deg', RA_deg, DE_deg));
            
            try

                total_time = 1;
                time_res = 0.1;
                N = ceil(total_time./time_res);
                total_time = N.*time_res; % if the division was not integer

                if obj.prechecks==0
                    obj.stop;
                    return;
                end

                obj.brake_bit = 0;
                on_cleanup = onCleanup(@obj.stop);
            
                % convert the total adjustment to the rate in arcsec/sec
                obj.hndl.RightAscensionRate = DE_deg.*3600/total_time;
                obj.hndl.DeclinationRate = DE_deg.*3600/total_time;
            
                t = tic;
                
                for ii = 1:N
                
                    if obj.brake_bit || obj.checks==0
                        obj.emergency_stop;
                        return;
                    end
                    
                    if toc(t)>total_time
                        obj.stop;
                        return;
                    end
                    
                    pause(time_res);
                    
                end
            
            catch ME
                obj.stop;
                obj.log.error(ME.getReport);
                rethrow(ME);
            end
                
        end
        
        function engineeringSlew(obj,Alt,Az)
            
            try 

                if Alt<obj.limit_alt
                    error('Input altitude is below limit of %f degress', obj.limit_alt);
                end
                
                obj.brake_bit = 0;
                
                obj.tracking = 0;
                obj.hndl.SlewToAltAzAsync(Alt, Az);
                
                for ii = 1:100000
                    
                    pause(0.01);
                    
                    if obj.brake_bit
                        obj.stop;
                        break;
                    end
                    
                    if ~obj.check_while_moving
                        obj.emergency_stop;
                        break;
                    end
                    
                    if obj.hndl.Slewing==0
                        break;
                    end
                    
                end
            
            catch ME
                obj.log.error(ME.getReport);
                rethrow(ME);
            end
            
            
        end
        
        function park(obj)
            
        end
        
        function zenith_west(obj)
            
        end
        
        function zenith_east(obj)
            
        end
        
        function update(obj)
            
            obj.target.update;
            
            try
                obj.hndl.Connected;
            catch 
                obj.status = 0; 
                return;
            end
            
            % check server is running
            [~, str] = system('tasklist /fi "imagename eq AstroOptikServer.exe"'); % check 
                
            str = strip(str);
            idx = strfind(str, 'AstroOptikServer.exe');
            
            if isempty(idx)
                obj.status = 0;
                return;
            end
            
            if obj.ALT<-20 % this occurs when software is disconnected from mount (e.g., on power out)
                obj.status = 0;
                return;
            end
            
            try
               
                if isempty(obj.ard) && obj.use_accelerometer
                    obj.ard = obs.sens.ScopeAssistant;
                    obj.ard.telescope = obj;
                end
                
                if ~isempty(obj.ard) && obj.use_accelerometer
                    obj.ard.update;
                end
                
            catch ME
                warning(ME.getReport);
            end
            
            % add additional tests?
            
            obj.status = 1;
            
        end
        
        function stop(obj)
            
%             obj.log.input('stopping telescope');
            
            try 
                
                obj.brake_bit = 1;
            
                obj.hndl.AbortSlew;
                
            catch ME
                obj.log.error(ME.getReport);
                rethrow(ME);
            end
            
        end
        
        function emergency_stop(obj)
            
            obj.log.input('stopping telescope');
            
            try 
                
                obj.brake_bit = 1;
            
                obj.hndl.AbortSlew;
                
            catch ME
                obj.log.error(ME.getReport);
                rethrow(ME);
            end
            
        end
        
    end
    
    methods % plotting tools / GUI
        
    end    
    
    methods(Static=true)
        
        
        
    end
    
end

