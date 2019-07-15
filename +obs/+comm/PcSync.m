classdef PcSync < handle


    properties(Transient=true)
        
        hndl_rx;
        hndl_tx;
        
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
        remote_port_rx = 4012;
        remote_port_tx = 4013;
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
        
        default_server_remote_port_rx = 4012;
        default_server_remote_port_tx = 4013;
        default_client_remote_port_rx = 4013;
        default_client_remote_port_tx = 4012;
        
        version = 1.01;
        
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
                        obj.remote_port_rx = obj.default_server_remote_port_rx;
                        obj.remote_port_tx = obj.default_server_remote_port_tx;
                        
                    elseif cs(varargin{1}, 'client')
                        
                        obj.role = 'Client';
                        obj.remote_ip = obj.default_client_remote_ip;
                        obj.remote_port_rx = obj.default_client_remote_port_rx;
                        obj.remote_port_tx = obj.default_client_remote_port_tx;
                        
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
        
        function val = get.is_connected(obj)
            
            if isempty(obj.hndl_rx) || isempty(obj.hndl_tx) ...
                    || ~isvalid(obj.hndl_rx) || ~isvalid(obj.hndl_tx) ...
                    || ~strcmp(obj.hndl_rx.Status, 'open') || ~strcmp(obj.hndl_tx.Status, 'open')
                val = 0;
            else
                val = 1;
            end
            
        end
        
    end
    
    methods % setters
        
    end
    
    methods % calculations
        
        function connect(obj, timeout)
            
            if nargin<2 || isempty(timeout)
                timeout = 10; % we will need to be smart about implementing this timeout, maybe using an extra worker
            end
            
            obj.log.input(['Setting up new connection to ' obj.remote_ip ':' num2str(obj.remote_port_tx) '/' num2str(obj.remote_port_rx) ' with role: ' obj.role]);
            
            try
                
                if util.text.cs(obj.role, 'server')
                    fprintf('Connecting to TCP/IP as server. There is no timeout! To break out hit Ctrl+C\n');
                end
                
                obj.disconnect;
                
                obj.hndl_tx = tcpip(obj.remote_ip, obj.remote_port_tx, 'NetworkRole', obj.role, 'Timeout', 10);
                obj.hndl_tx.BytesAvailableFcn = @obj.read_data;
                obj.hndl_tx.BytesAvailableFcnMode = 'byte';
                obj.hndl_tx.BytesAvailableFcnCount = 32;
                obj.hndl_tx.OutputBufferSize = 50*1024; 
                obj.hndl_tx.InputBufferSize = 50*1024;
                
                fopen(obj.hndl_tx);
            
                obj.hndl_rx = tcpip(obj.remote_ip, obj.remote_port_rx, 'NetworkRole', obj.role, 'Timeout', 10);
                obj.hndl_rx.BytesAvailableFcn = @obj.read_data;
                obj.hndl_rx.BytesAvailableFcnMode = 'byte';
                obj.hndl_rx.BytesAvailableFcnCount = 32;
                obj.hndl_rx.OutputBufferSize = 50*1024; 
                obj.hndl_rx.InputBufferSize = 50*1024;
                
                fopen(obj.hndl_rx);
            
                obj.update;
                
            catch ME
                obj.log.error(ME.getReport);
                rethrow(ME);
            end
            
        end
        
        function disconnect(obj)
            
            if ~isempty(obj.hndl_tx) && isvalid(obj.hndl_tx)
                delete(obj.hndl_tx);
                obj.hndl_tx = [];
            end
            
            if ~isempty(obj.hndl_rx) && isvalid(obj.hndl_rx)
                delete(obj.hndl_rx);
                obj.hndl_rx = [];
            end

        end
        
        function send(obj, value, rx_or_tx)
            
            if nargin<3 || isempty(rx_or_tx)
                rx_or_tx = 'tx';
            end
            
            if strcmpi(rx_or_tx, 'tx') % primary transmission mode
                obj.raw_data_sent = getByteStreamFromArray(value);
                obj.checksum = util.oop.getHash(obj.raw_data_sent);
                fwrite(obj.hndl_tx, obj.raw_data_sent);
            elseif strcmpi(rx_or_tx, 'rx') % reply only (e.g., sending back the hash of latest incoming data)
                temp_raw_data = getByteStreamFromArray(value);
                fwrite(obj.hndl_rx, temp_raw_data);
            else
                error('Must choose RX or TX for 3rd input to send()');
            end
                
        end
        
        function update(obj)
            
            if obj.is_connected
                obj.hndl_rx.BytesAvailableFcn = @obj.read_data;
                obj.outgoing.time = util.text.time2str(datetime('now', 'TimeZone', 'UTC'));

                obj.status = 0;
                flushinput(obj.hndl_tx);
                obj.send(obj.outgoing);
            end
            
        end
        
        function reply_hash(obj)
            
            obj.send(util.oop.getHash(obj.raw_data_received), 'rx');
            
        end
        
        function read_data(obj, hndl, ~)
            
            if obj.debug_bit>1, fprintf('read data with %d bytes\n', hndl.BytesAvailable); end
            
            if hndl.BytesAvailable>0
                
                data = uint8(fread(hndl, hndl.BytesAvailable))';
                
                try
                    value = getArrayFromByteStream(data);
                catch ME
                    value = [];
                    warning(ME.getReport)
                end
                
                if isempty(value)
                    % pass 
                elseif ischar(value)
                    if strcmp(obj.checksum, value)
                        obj.status = 1;
                    else
                        error('Received a response: %s which is not consistent with checksum: %s', value, obj.checksum);
                    end
                elseif isstruct(value)
                    obj.raw_data_received = data;
                    obj.incoming = value;
                    obj.status = 1;
                    obj.reply_hash;
                end
                
            end
            
        end
        
    end
    
    methods % plotting tools / GUI
        
    end    
    
end

