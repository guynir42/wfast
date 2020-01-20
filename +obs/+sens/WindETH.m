% WindETH class  
% Package: +obs/+sens/
% Description: A class for controlling the WindETH ethernet anemometer
%              (from papouch.com).
%              https://www.papouch.com/en/shop/product/windeth-ethernet-anemometer/
%              The class open a tcp/ip port to the instrument
% Tested : Matlab R2018a
%     By : Eran O. Ofek                    Nov 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: % Create a WindETH object
%          WC=obs.sens.WindETH('132.77.38.187',10001);
%          % update wind speed/az values in WC object properties
%          getWind(WC);
%          % close/open WC.ComObj tcpip object
%          open(WC); close(WC);
% Reliable: 2
%--------------------------------------------------------------------------

classdef WindETH < handle
    
    properties % objects
        
        hndl;
        
    end
    
    properties
        
        % generic fields
        id = 'wind';
        status = 0; % false - readings are unreliable, true - ok
        use_this = 1; % is true get the data from this sensor. If not, ignore it and move on
        
        data@util.vec.CircularBuffer; % Stack object containing data history
        dataCol = {'JD','WindSpeed','WindAz'};   % Stack object columns
        dataUnits = {'day','km/h','deg'};          % Stack object column units
        
        ip = '192.168.1.254';
        port = 10001;
        
        % replies from hndl
        jd = NaN; % JD of last sucessful reading
        wind_speed = NaN; % last wind speed
        wind_speed_status; % if speed measurement is ok or error
        wind_az = NaN; % last wind azimuth
        wind_az_status; % if azimuth measurement is ok or error
        
        reply = ''; % verbatim reply from hardware
        
    end
    
    properties (Hidden = true)
        
        north_az = 112.5;
        
        version = 1.01;
        
    end
    
    methods % Constructor
        
        function obj=WindETH(ip,port)
        % Example: obj = obs.sens.WindETH('132.77.38.187',10001)
            
            if nargin>=1 && ~isempty(ip)
                obj.ip = ip;
            end
            
            if nargin>=2 && ~isempty(port)
                obj.port = port;
            end
            
            obj.data = util.vec.CircularBuffer(100);
            
            obj.connect(obj.ip, obj.port);
            
            obj.update;
            
        end
        
        function connect(obj, ip, port)
            
            if nargin>1 && ~isempty(ip)
                obj.ip = ip;
            end
            
            if nargin>2 && ~isempty(port)
                obj.port = port;
            end
            
            obj.disconnect;
            
            obj.hndl = tcpip(obj.ip, obj.port);
            
            % obj.hndl.Timeout = 10;
            obj.hndl.Terminator = 'CR';  % set terminaotor to carige return
            
            fopen(obj.hndl);
                        
            
        end
        
        function disconnect(obj)
            
            if ~isempty(obj.hndl) && isvalid(obj.hndl)
                fclose(obj.hndl);
                delete(obj.hndl);
                obj.hndl = [];
            end
            
        end
        

    end
    
    methods % resetters
        
        function reset(obj)
            
            obj.data.reset;
            
        end
        
    end
    
    methods 
        
        function update(obj)
            
            obj.status = 0;
            
            % send query to WindETH and get answer
            %fwrite(obj.hndl,unicode2native(sprintf('*B1MR0\r'), 'UTF-8'))
            %V=fread(obj.hndl);
            % for query see Papouch s.r.o Spinel in AD4xxx manual: pp 10-11
            obj.reply = query(obj.hndl, sprintf('*B1MR0\r'));
            
            if isempty(obj.reply)
                error('Cannot get a reply from WindETH');
            end
                        
            % single measurment command
%             HexCommand = ['2A';'61';'00';'06';'31';'02';'51';'00';'EA';'0D'];
%             DecCommand = hex2dec(HexCommand); 
%             fwrite(WC.ComObj,DecCommand,'uint8');
%             DecAns     = fread(WC.ComObj);
%             HexAns     = dec2hex(DecAns); 
%             
%             HexAns(7,:)     % command recieved OK
%             HexAns(8,:)     % status for ch1
%             HexAns(9:10,:)  % ch1 value
%             HexAns(12,:)    % status for ch2
%             HexAns(13:14,:) % ch2 value
%             
            % Extract data from output string 
            Cell = regexp(obj.reply,'\s','split');
            Flag = cellfun(@isempty,Cell);
            Cell = Cell(~Flag);
            
            WindAzStatus    = Cell{3};   % WindAz status flag [80 | 88]
            WindAzString    = Cell{4};   % WindAz direction string
            WindSpeedStatus = Cell{6};   % WindSpeed status flag [80 | 88]
            WindSpeedString = Cell{7};   % WindSpeed 
            
            obj.status = 1;
            
            % set Wind speed/dir status
            switch WindSpeedStatus
                case '80'
                    obj.wind_speed_status = 'ok';
                case '88'
                    obj.wind_speed_status = 'out of range';
                otherwise
                    obj.wind_speed_status = 'error';
                    obj.status = 0;
            end
            
            switch WindAzStatus
                case '80'
                    obj.wind_az_status = 'ok';
                case '88'
                    obj.wind_az_status = 'out of range';
                otherwise
                    obj.wind_az_status = 'error';
                    obj.status = 0;
            end
            
            % convert speed/dir to numeric values
            if (obj.status)
                % if status ok - convert wind Az/spped to numeric values
                obj.wind_speed = str2double(WindSpeedString);
                
                Az = obs.sens.WindETH.dir2az(WindAzString);
                
                if ~isempty(obj.north_az)
                    Az = mod(Az - obj.north_az, 360);
                end
                
                if (isnan(Az))
                    obj.status = 0;
                end
                
                obj.wind_az = Az;
                
            end
            
            % if status is ok than set the time of the last query
            if (obj.status)
               % update time to that of last measurment
               obj.jd = juliandate(datetime('now', 'timezone', 'utc'));
            end
            
            obj.data.input([obj.jd, obj.wind_speed, obj.wind_az]);
             
        end
        
    end
    
    % static methods
    methods (Static)
        
        function Az=dir2az(DirStr)
            % Convert direction (e.g., NEE) to Azimuth
            % Input  : Direction (e.g., N, NE, NNW)
            % Output : Azimuth (deg)
            % Example: obs.sens.WindETH.dir2az('NE')
            
            N_N = numel(strfind(DirStr,'N'));
            N_E = numel(strfind(DirStr,'E'));
            N_S = numel(strfind(DirStr,'S'));
            N_W = numel(strfind(DirStr,'W'));
            
            N = [N_N N_E N_S N_W];
            
            switch sprintf('%d%d%d%d',N)
                case '1000'
                    Az = 0;
                case '0100'
                    Az = 90;
                case '0010'
                    Az = 180;
                case '0001'
                    Az = 270;
                case '2100'
                    Az = 22.5;
                case '1100'
                    Az = 45;
                case '1200'
                    Az = 67.5;
                case '0210'
                    Az = 112.5;
                case '0110'
                    Az = 135;
                case '0120'
                    Az = 157.5;
                case '0021'
                    Az = 202.5;
                case '0011'
                    Az = 225;
                case '0012'
                    Az = 247.5;
                case '1002'
                    Az = 292.5;
                case '1001'
                    Az = 315;
                case '2001'
                    Az = 337.5;
                otherwise
                    Az = NaN;
             end
            
        end
        
    end
    
end

            
