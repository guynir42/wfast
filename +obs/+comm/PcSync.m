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
        
        remote_ip = '192.168.1.103';
        remote_port_rx = 4012;
        remote_port_tx = 4013;
        role = '';
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        is_connected;
        
    end
    
    properties(Hidden=true)
       
        raw_data_rx_temp; % for broken up incoming messages
        raw_data_tx_temp;
        
        default_client_remote_ip = '192.168.1.101';
        default_server_remote_ip = '0.0.0.0';
        
        default_server_remote_port_rx = 4012;
        default_server_remote_port_tx = 4013;
        default_client_remote_port_rx = 4013;
        default_client_remote_port_tx = 4012;
        
        version = 1.02;
        
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
                
                if strcmpi(obj.role, 'server')
                    obj.hndl_rx = obj.connectSocket('rx');
                    obj.hndl_tx = obj.connectSocket('tx');
                else
                    obj.hndl_tx = obj.connectSocket('tx');
                    obj.hndl_rx = obj.connectSocket('rx');
                end
                
                pause(0.1);
                
                obj.update;
                
            catch ME
                obj.log.error(ME.getReport);
                rethrow(ME);
            end
            
        end
        
        function hndl = connectSocket(obj, rx_or_tx)

            hndl = tcpip(obj.remote_ip, obj.(['remote_port_' rx_or_tx]), 'NetworkRole', obj.role, 'Timeout', 10);
            hndl.OutputBufferSize = 50*1024; 
            hndl.InputBufferSize = 50*1024;
            flushinput(hndl);
            flushoutput(hndl);

            try
                fopen(hndl);
            catch
                delete(hndl); 
            end
            
%             hndl.BytesAvailableFcnCount = 32;
            hndl.BytesAvailableFcnMode = 'terminator';
            if strcmp(rx_or_tx, 'tx')
                hndl.BytesAvailableFcn = @obj.read_data_tx;
            else
                hndl.BytesAvailableFcn = @obj.read_data_rx;
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
            
            byte_stream = getByteStreamFromArray(value);
            
            if strcmpi(rx_or_tx, 'tx') % primary transmission mode
                obj.raw_data_sent = uint8(sprintf('%s Message length is %010d\n', byte_stream, length(byte_stream)));
                obj.checksum = util.oop.getHash(obj.raw_data_sent);
                obj.waitForTransferStatus(obj.hndl_tx);
                fwrite(obj.hndl_tx, obj.raw_data_sent);
            elseif strcmpi(rx_or_tx, 'rx') % reply only (e.g., sending back the hash of latest incoming data)
                obj.waitForTransferStatus(obj.hndl_rx);
                fwrite(obj.hndl_rx, sprintf('%s Message length is %010d\n', byte_stream, length(byte_stream))); 
            else
                error('Must choose RX or TX for 3rd input to send()');
            end
                
        end
        
        function waitForTransferStatus(~, hndl)
            
            res = 0.01;
            
            timeout = 5;
            
            for ii = 1:timeout/res
                
                if strcmp(hndl.TransferStatus, 'idle')
                    return;
                end
                
                pause(res);
                
            end
            
            error('Timeout afer %f seconds while waiting for tcp/ip object to clear the TransferStatus', timeout);
            
        end
        
        function update(obj)
            
            if obj.is_connected
                
                obj.hndl_rx.BytesAvailableFcn = @obj.read_data_rx;
                obj.hndl_tx.BytesAvailableFcn = @obj.read_data_tx;
                obj.outgoing.time = util.text.time2str(datetime('now', 'TimeZone', 'UTC'));

                obj.status = 0;
%                 flushinput(obj.hndl_tx);
                obj.send(obj.outgoing);
                
            end
            
        end
        
        function reply_hash(obj)
            
            obj.send(util.oop.getHash(obj.raw_data_received), 'rx');
            
        end
        
        function read_data(obj, hndl, rx_or_tx) % can be called separately for rx or tx
            
            if obj.debug_bit>1, fprintf('Received data with %d bytes\n', hndl.BytesAvailable); end
            
%             obj.waitForTransferStatus(hndl);
            
            variable = [];
            
            data_name = ['raw_data_' rx_or_tx '_temp'];
            
            data_temp = obj.(data_name);
            
            for ii = 1:10

                if hndl.BytesAvailable>0
                    
                    obj.waitForTransferStatus(hndl);
                    
                    new_data = uint8(fread(hndl, hndl.BytesAvailable))'; % this may be only a part of the message...
                                       
                    data_temp = horzcat(data_temp, new_data); % add the latest data to temp buffer
                    
                    if length(data_temp)<30 % message is too short (no footer attached!)
                        if obj.debug_bit>1, disp(['Message is too short, only ' num2str(length(data_temp))]); end
                        break;
                    end
                    
                    footer_str = char(data_temp(end-28:end-1)); % include only the text and number (no line break, no leading space)
                    
                    if ~strcmp(footer_str(1:18), 'Message length is ') % message doesn't have correct footer string
                        if obj.debug_bit>1, disp(['Incorrect footer string: ' footer_str]); end
                        break;
                    end
                                        
                    idx = regexp(char(footer_str), '\d{10}');
                    
                    if isempty(idx) || idx~=19
                        if obj.debug_bit>1, disp(['Numeric values index is ' num2str(idx) '... should be at 19']); end
                        break;
                    end
                    
                    L = str2double(footer_str(19:end));
                    
                    if length(data_temp)<L
                        if obj.debug_bit>1, fprintf('Length of message is %d, smaller than required by footer data (%d bytes). Dropping entire message!\n', length(data_temp), L); end
                        obj.(data_name) = uint8([]);
                        return;
                    end
                    
                    try
                    
                        variable = getArrayFromByteStream(data_temp(end-30-L+1:end-30)); % read L bytes from collected data, minus footer string. 
                        break; % upon success, should break from the loop also... 
                    
                    catch ME % failed to deserialize the values... 
                        
                        if strcmp(ME.identifier, 'MATLAB:Deserialize:BadVersionOrEndian')
                            
                            if obj.debug_bit>1
                                disp(['"data" cannot be parsed after ' num2str(length(data_temp)-30) ' bytes... try to append more!']); % this should not happen any more! 
                            end
                            
                            continue;
                            
                        else
                            disp(['Got a different error: ' ME.identifier]);
                            warning(ME.getReport);
                            obj.(data_name) = uint8([]); % upon a generic error, should drop the whole message
                            return;
                        end
                        
                    end
                    
                end
                
            end % for ii 
            
            obj.(data_name) = data_temp; % if no "return" was issued, update the raw_data_rx_temp or raw_data_tx_temp with latest data
            
            if isempty(variable)
                if obj.debug_bit>1, disp('Received an empty variable'); end
            elseif strcmp(rx_or_tx, 'tx') && ischar(variable)
                if ~strcmp(obj.checksum, variable)
                    fprintf('Time: %s | message length= %d | pipe: %s\n', util.text.time2str(datetime('now', 'TimeZone', 'UTC')), length(data_temp), rx_or_tx); 
                    warning('Checksum for latest transmission was %s, received confirmation checksum: %s', obj.checksum, variable)
                end
                obj.status = 1;
            elseif strcmp(rx_or_tx, 'rx') && isstruct(variable)
                if obj.debug_bit>1, disp(['Successfully deserialized message with ' num2str(length(obj.raw_data_rx_temp)-1) ' bytes! Converted into a struct... ']); end
                obj.raw_data_received = data_temp;
                obj.(data_name) = uint8([]);
                obj.incoming = variable;
                obj.status = 1;
                obj.reply_hash;
            else
                warning('Wrong type of data received in %s pipe, class(variable)= %s', rx_or_tx, class(variable)); % this shouldn't happen! 
            end

        end
        
        function read_data_rx(obj, hndl, ~)
            
            obj.read_data(hndl, 'rx');

        end
        
        function read_data_tx(obj, hndl, ~)
            
            obj.read_data(hndl, 'tx');

        end
        
    end
    
    methods % plotting tools / GUI
        
    end    
    
end

