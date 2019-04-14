classdef AudioControl <handle
% class that plays sound effects (e.g. when starting/ending analysis)
%
% TEST PROTOCOL: a = util.sys.AudioControl; a.playTakeForever;
    
    properties(SetAccess='protected')
       
        % sound effects
        horn;
        siren;
        dubstep;
        shows_over;
        take_forever;
        alarm_alarm;
        black_alert;
        rick_morty;
        what_you_got;
        warning;
        % add the other sound effects...
        
    end
    
    properties
        
        master_switch = 1;
        
        max_duration = 25; % seconds
        
        stop_timer = timer;
        
        debug_bit = 1;
        
    end
    
    properties(Hidden=true)
       
        version = 1.02;
        
    end
    
    methods % constructor
       
        function obj = AudioControl
           
            if obj.debug_bit, fprintf('AudioControl constructor v%4.2f\n', obj.version); end
            
            try 
                
                [Y, Fs] = audioread([getenv('WFAST') '/+util/audio_clips/horn_regular.mp3']);
                obj.horn = audioplayer(Y, Fs);

                [Y, Fs] = audioread([getenv('WFAST') '/+util/audio_clips/horn_siren.mp3']);
                obj.siren = audioplayer(Y, Fs);
            
                [Y, Fs] = audioread([getenv('WFAST') '/+util/audio_clips/alarm_dubstep.mp3']);
                obj.dubstep = audioplayer(Y, Fs);
            
                [Y, Fs] = audioread([getenv('WFAST') '/+util/audio_clips/takeForever.mp3']);
                obj.take_forever = audioplayer(Y, Fs);

                [Y, Fs] = audioread([getenv('WFAST') '/+util/audio_clips/showsOver.mp3']);
                obj.shows_over = audioplayer(Y, Fs);

                [Y, Fs] = audioread([getenv('WFAST') '/+util/audio_clips/alarm_alarm.mp3']);
                obj.alarm_alarm = audioplayer(Y, Fs);

                [Y, Fs] = audioread([getenv('WFAST') '/+util/audio_clips/black_alert.mp3']);
                obj.black_alert = audioplayer(repmat(Y*5, [2,1]), Fs);
                
                [Y, Fs] = audioread([getenv('WFAST') '/+util/audio_clips/rick_morty.mp3']);
                obj.rick_morty = audioplayer(Y, Fs);
                
                [Y, Fs] = audioread([getenv('WFAST') '/+util/audio_clips/what_you_got.mp3']);
                obj.what_you_got = audioplayer(Y*3, Fs);
                
                [Y, Fs] = audioread([getenv('WFAST') '/+util/audio_clips/warning.mp3']);
                obj.warning = audioplayer(Y*5, Fs);
                
            catch ME
                warning(ME.getReport);
                obj.master_switch = 0;
            end
            
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
            
            obj.stop_timer = timer;
            
            set(obj.stop_timer, 'Name','audio-stop-timer');
%             set(obj.stop_timer, 'ExecutionMode','fixedRate');
            set(obj.stop_timer, 'StartDelay',obj.max_duration);
            set(obj.stop_timer, 'TimerFcn',{@obj.stop});
        
            start(obj.stop_timer);
            
        end
                
        function stopTimer(obj)
            
            stop(obj.stop_timer);
            delete(obj.stop_timer);
            
        end
            
    end
    
    methods % audio control
        
        function stop(obj, ~, ~)
            
            try
            
                obj.horn.pause;
                obj.siren.pause;
                obj.dubstep.pause;
                obj.shows_over.pause;
                obj.take_forever.pause;
                obj.alarm_alarm.pause;
                obj.black_alert.pause;
                obj.rick_morty.pause;
                obj.what_you_got.pause;
                obj.warning.pause;
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
                
                obj.setupTimer
            
                obj.horn.play;
                
            catch ME
                warning(ME.getReport);
            end
            
        end
        
        function playSiren(obj)
        
            try
                
                if obj.master_switch==0
                    return;
                end
                
                obj.setupTimer
            
                obj.siren.play;
                
            catch ME
                warning(ME.getReport);
            end
            
        end
        
        function playDubstep(obj)
        
            try
                
                if obj.master_switch==0
                    return;
                end
                
                obj.setupTimer
            
                obj.dubstep.play;
                
            catch ME
                warning(ME.getReport);
            end
            
        end
        
        function playShowsOver(obj)
            
            try

                if obj.master_switch==0
                    return;
                end
                
                obj.setupTimer
            
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
                
                obj.setupTimer
            
                obj.take_forever.play;
            
            catch ME
                warning(ME.getReport);
            end
            
        end
        
        function playAlarmAlarm(obj)
            
            try
            
                if obj.master_switch==0
                    return;
                end
                
                obj.setupTimer
            
                obj.alarm_alarm.play;
            
            catch ME
                warning(ME.getReport);
            end
            
        end
        
        function playBlackAlert(obj)
            
            try
            
                if obj.master_switch==0
                    return;
                end
                
                obj.setupTimer
            
                obj.black_alert.play;
            
            catch ME
                warning(ME.getReport);
            end
            
        end
        
        function playRickMorty(obj)
            
            try
            
                if obj.master_switch==0
                    return;
                end
                
                obj.setupTimer
            
                obj.rick_morty.play;
            
            catch ME
                warning(ME.getReport);
            end
            
        end
        
        function playWhatYouGot(obj)
            
            try
            
                if obj.master_switch==0
                    return;
                end
                
                obj.setupTimer
            
                obj.what_you_got.play;
            
            catch ME
                warning(ME.getReport);
            end
            
        end
        
        function playWarning(obj)
            
            try
            
                if obj.master_switch==0
                    return;
                end
                
                obj.setupTimer
            
                obj.warning.play;
            
            catch ME
                warning(ME.getReport);
            end
            
        end
        
        function playError(obj)
           
            % not yet implemented
            
        end
        
    end
    
end