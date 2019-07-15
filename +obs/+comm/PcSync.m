classdef PcSync < handle


    properties(Transient=true)
        
        hndl;
        
    end
    
    properties % objects
        
        incoming; % struct with data that got received from other side
        outgoing; % struct with data to be sent to the other side
        
        log@util.sys.Logger;
        
    end
    
    properties % inputs/outputs
        
        raw_data_received;
        raw_data_sent;
        checksum; % for sent data
        
        status = 0;
        
    end
    
    properties % switches/controls
        
        remote_ip = '192.168.1.100';
        remote_port = 4012;
        role = '';
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        name;
        timeout;
        is_connected;
        
    end
    
    properties(Hidden=true)
       
        default_client_remote_ip = '192.168.1.101';
        default_server_remote_ip = '192.168.1.100';
        
        version = 1.00;
        
    end
    
    methods % constructor
        
        function obj = PcSync(varargin)
            
            import util.text.cs;
            
            if ~isempty(varargin) && isa(varargin{1}, 'obs.comm.PcSync')
                if obj.debug_bit, fprintf('PcSync copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit, fprintf('PcSync constructor v%4.2f\n', obj.version); end
            
                if isscalar(varargin) && ischar(varargin{1})
                    
                    if cs(varargin{1}, 'server')
                        obj.role = 'Server';
                        obj.remote_ip = obj.default_server_remote_ip;
                    elseif cs(varargin{1}, 'client')
                        obj.role = 'Client';
                        obj.remote_ip = obj.default_client_remote_ip;
                    % add additional initializations? 
                    end
                    
                    if isempty(obj.role)
                        obj.log = util.sys.Logger('comm_sync');
                    else
                        obj.log = util.sys.Logger(['comm_sync_' obj.role]);
                    end
                    
                    obj.reset;
                    
                end
                
            end
            
        end
        
        function delete(obj)
            
            obj.disconnect;
            
        end

    end
    
    methods % reset/clear
        
        function reset(obj)
            
            obj.incoming = struct;
            obj.outgoing = struct;
            
        end
        
    end
    
    methods % getters
        
        function val = get.name(obj)
            
            if isempty(obj.hndl)
                val = [];
            else
                val = obj.hndl.Name;
            end
            
        end
        
        function val = get.timeout(obj)
            
            if isempty(obj.hndl)
                val = [];
            else
                val = obj.hndl.Timeout;
            end
            
        end
        
        function val = get.is_connected(obj)
            
            if isempty(obj.hndl) || ~isvalid(obj.hndl) || ~strcmp(obj.hndl.Status, 'open')
                val = 0;
            else
                val = 1;
            end
            
        end
        
    end
    
    methods % setters
        
        function set.name(obj, val)
            
            if ~isempty(obj.hndl)
                obj.hndl.Name = val;
            end
            
        end
        
        function set.timeout(obj, val)
            
            if ~isempty(obj.hndl)
                obj.hndl.Timeout = val;
            end
            
        end
        
    end
    
    methods % calculations
        
        function connect(obj, timeout)
            
            if nargin<2 || isempty(timeout)
                timeout = 10; % we will need to be smart about implementing this timeout, maybe using an extra worker
            end
            
            obj.log.input(['Setting up new connection to ' obj.remote_ip ':' num2str(obj.remote_port) ' with role: ' obj.role]);
            
            try
                
                if util.text.cs(obj.role, 'server')
                    fprintf('Connecting to TCP/IP as server. There is no timeout! To break out hit Ctrl+C\n');
                end
                
                obj.disconnect;
                
                obj.hndl = tcpip(obj.remote_ip, obj.remote_port, 'NetworkRole', obj.role, 'Timeout', 10);
                obj.hndl.BytesAvailableFcn = @obj.read_data;
                obj.hndl.BytesAvailableFcnMode = 'byte';
                obj.hndl.BytesAvailableFcnCount = 32;
                obj.hndl.OutputBufferSize = 50*1024; 
                obj.hndl.InputBufferSize = 50*1024;
                
                fopen(obj.hndl);
            
                obj.update;
                
            catch ME
                obj.log.error(ME.getReport);
                rethrow(ME);
            end
            
        end
        
        function disconnect(obj)
            
            if ~isempty(obj.hndl) && isvalid(obj.hndl)
                delete(obj.hndl);
                obj.hndl = [];
            end

        end
        
        function send(obj, value)
            
            obj.raw_data_sent = getByteStreamFromArray(value);
            obj.checksum = util.oop.getHash(obj.raw_data_sent);
            fwrite(obj.hndl, obj.raw_data_sent);
            
        end
        
        function update(obj)
            
            if obj.is_connected
                obj.hndl.BytesAvailableFcn = @obj.read_data;
                obj.outgoing.time = util.text.time2str(datetime('now', 'TimeZone', 'UTC'));

                obj.status = 0;
                flushinput(obj.hndl);
                obj.send(obj.outgoing);
            end
            
        end
        
        function confirm(obj)
            
            obj.send(util.oop.getHash(obj.raw_data_received));
            
        end
        
        function read_data(obj, ~, ~)
            
            if obj.debug_bit>1, fprintf('read data with %d bytes\n', obj.hndl.BytesAvailable); end
            
            if obj.hndl.BytesAvailable>0
                
                obj.raw_data_received = uint8(fread(obj.hndl, obj.hndl.BytesAvailable))';
                
                try
                    value = getArrayFromByteStream(obj.raw_data_received);
                catch ME
                    value = [];
                    warning(ME.getReport)
                end
                
                if isempty(value)
                    % what to do here??
                elseif ischar(value)
                    if strcmp(obj.checksum, value)
                        obj.status = 1;
                    else
                        error('Received a response: %s which is not consistent with checksum: %s', value, obj.checksum);
                    end
                elseif isstruct(value)
                    obj.incoming = value;
                    obj.status = 1;
                    obj.confirm;
                end
                
            end
            
        end
        
    end
    
    methods % plotting tools / GUI
        
    end    
    
end

