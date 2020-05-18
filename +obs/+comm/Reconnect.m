classdef Reconnect < handle

    properties(Transient=true)
        
    end
    
    properties % objects
        
    end
    
    properties % inputs/outputs
        
        state_of_lock = false; 
        start_time_lock; 
        num_attempts = 0;
        history = {};
        
    end
    
    properties % switches/controls
        
        delay_time_minutes = 30; % how many minutes to wait before unlocking
        max_failed_attempts = 2; % how many failed connections would lock
        
        size_history = 30;
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        
        
    end
    
    properties(Hidden=true)
       
        version = 1.00;
        
    end
    
    methods % constructor
        
        function obj = Reconnect(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'obs.comm.Reconnect')
                if obj.debug_bit>1, fprintf('Reconnect copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit>1, fprintf('Reconnect constructor v%4.2f\n', obj.version); end
            
            end
            
        end
        
    end
    
    methods % reset/clear
        
        function reset(obj)
            
            obj.unlock;
            obj.history = {};
            
        end
        
    end
    
    methods % getters
        
    end
    
    methods % setters
        
    end
    
    methods % calculations
        
        function val = should(obj) % check if it OK to try again to reconnect to device
            
            if obj.state_of_lock
                
                % test if enough time has passed
                if etime(clock, obj.start_time_lock)/60 > obj.delay_time_minutes
                    obj.unlock;
                    val = true; % enough time has passed to unlock and try again
                else
                    val = false; % we are still locked! 
                end
                
            else
                val = true; % if unlocked, just let them reconnect
            end
            
        end
        
        function lock(obj) % too many attempts have been made, locking this object for a few minutes
            
            obj.start_time_lock = clock;
            obj.state_of_lock = true;
            
        end
        
        function unlock(obj) % no longer need to be locked
        
            obj.state_of_lock = false;
            obj.num_attempts = 0;
            
        end
        
        function inputSuccess(obj) % successfully connected to device, can reset all failed attempts
            
            obj.unlock;
            obj.addToHistory('Successful reconnect!'); 
            
        end
        
        function inputFailure(obj, report) % failed to connect to device, maybe time to lock this object? optional argument added to history
            
            if nargin<2 || isempty(report)
                report = '';
            end
            
            obj.num_attempts = obj.num_attempts + 1;
            
            if obj.num_attempts > obj.max_failed_attempts
                obj.lock;
                obj.addToHistory(sprintf('Failed to connect to device %d times! Locking for %d minutes...', obj.num_attempts, obj.delay_time_minutes));
            else
                obj.addToHistory(sprintf('Failed to connect to device %d times!', obj.num_attempts)); 
            end
            
            if ~isempty(report)
%                 obj.addToHistory(util.text.eraseTags(report)); 
                obj.addToHistory(report); 
            end
            
        end
        
        function addToHistory(obj, str) % add a new string to the history cell array (with date pre-appended) 
            
            t = util.text.time2str(datetime('now', 'TimeZone', 'UTC')); 
            
            obj.history{end+1} = sprintf('%s: %s', t, str);
            
            if size(obj.history,2)>1
                obj.history = obj.history';
            end
            
            if length(obj.history)>obj.size_history
                obj.history = obj.history(end-obj.size_history+1:end); 
            end
            
        end
            
    end
    
    methods % plotting tools / GUI
        
    end    
    
end

