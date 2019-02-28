classdef AudioControl <handle
% class that plays sound effects (e.g. when starting/ending analysis)
%
% TEST PROTOCOL: a = util.sys.AudioControl; a.playTakeForever;
    
    properties(SetAccess='protected')
       
        % sound effects
        horn;
        shows_over;
        take_forever;
        % add the other sound effects...
        
    end
    
    properties
        
        master_switch = 1;
        
        max_duration=5; % seconds
        
        stop_timer = timer;
        
        debug_bit = 1;
        
    end
    
    properties(Hidden=true)
       
        version = 1.01;
        
    end
    
    methods % constructor
       
        function obj = AudioControl
           
            if obj.debug_bit, fprintf('AudioControl constructor v%4.2f\n', obj.version); end
            
            try 
            [Y, Fs] = audioread([getenv('MAT') '/+util/audio_clips/horn.mp3']);
            obj.horn = audioplayer(Y, Fs);
            
            [Y, Fs] = audioread([getenv('MAT') '/+util/audio_clips/takeForever.mp3']);
            obj.take_forever = audioplayer(Y, Fs);
            
            [Y, Fs] = audioread([getenv('MAT') '/+util/audio_clips/showsOver.mp3']);
            obj.shows_over = audioplayer(Y, Fs);
            
            catch ME
                warning(ME.getReport);
                master_switch = 0;
            end
            
            obj.setupTimer;
            
        end
        
    end
    
    methods % setters
       
        function set.max_duration(obj, duration)
           
            if duration<0
                obj.max_duration = 0;
            else
                obj.max_duration = duration;
            end
            
        end
        
    end
        
    methods % timer methods
        
        function setupTimer(obj)
            
            stop(obj.stop_timer);
            
            set(obj.stop_timer, 'Name','audio-stop-timer');
            set(obj.stop_timer, 'ExecutionMode','fixedRate');
            set(obj.stop_timer, 'Period',obj.max_duration);
            set(obj.stop_timer, 'TimerFcn',{@obj.stop});
            
            obj.startTimer;
            
        end
        
        function startTimer(obj)
            
            if util.text.cs(obj.stop_timer.Running, 'on')
                stop(obj.stop_timer);
            end
            
            start(obj.stop_timer);
            
        end
        
        function stopTimer(obj)
            
            stop(obj.stop_timer);
            
        end
            
    end
    
    methods % audio control
        
        function stop(obj, hndl, event)
           
            if nargin==1
                
            end
            
            try
            
                obj.horn.pause;
                obj.take_forever.pause;
                obj.shows_over.pause;
                % add the other sound effects...

                obj.stopTimer;

            catch ME
                warning(ME.getReport);
            end
            
        end
        
        function playHorn(obj)
        
            try
                
                if obj.master_switch==0
                    return;
                end
                
                obj.startTimer;
            
                obj.horn.play;
                
            catch ME
                warning(ME.getReport);
            end
            
        end
        
        function playShowsOver(obj)
            
            try

                if obj.master_switch==0
                    return;
                end
                
                obj.startTimer;
            
                obj.shows_over.play;
                           
            catch ME
                warning(ME.getReport);
            end
            
        end
        
        function playTakeForever(obj)
            
            try
            
                if obj.master_switch==0
                    return;
                end
                
                obj.startTimer;
            
                obj.take_forever.play;
            
            catch ME
                warning(ME.getReport);
            end
            
        end
        
        function playError(obj)
           
            % not yet implemented
            
        end
        
    end
    
end