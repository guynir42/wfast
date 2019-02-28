classdef ProgressBar < handle
% Keeps track of time and prints out progress. 
%
% TEST PROTOCOL: p=util.ProgressBar(100); p.showif(10); p.showif(100);

    properties
        
        start_time;
        stored_time;
        is_paused = 1;
        
        total_number;
        current_number;
        
        print_bit = 1;
                        
    end
    
    properties(Hidden=true)
        
    end
    
    properties(Dependent=true)
        
        dividor;
        
    end
    
    methods % constructor
        
        function obj = ProgressBar(total_number)
                                    
            if nargin>0 && ~isempty(total_number)
                obj.total_number = total_number;
            end
            
            obj.reset;
                     
            
        end
        
    end
    
    methods % getters
       
        function d = get.dividor(obj)
           
            d = 10^floor(log10(obj.total_number)-0.5);
            
%             if d<10
%                 d = 10;
%             end
            
        end
        
        function t = getElapsed(obj)
           
            if obj.is_paused
                t = obj.stored_time;
            else
                t = etime(clock, obj.start_time) + obj.stored_time;
            end
            
        end
        
    end
    
    methods % actions
                
        function start(obj, total_number)
               
            if nargin>1 && ~isempty(total_number)
                obj.total_number = total_number;
            end
            
            obj.reset(total_number);
            obj.unpause;
            
        end
        
        function finish(obj, number)
        
            if nargin<2
                number = [];
            end
            
            obj.show(number);
            obj.pause;
            
        end
        
        function reset(obj, total_number)
            
            if nargin>1 && ~isempty(total_number)
                obj.total_number = total_number;
            end

            obj.stored_time = 0;
            obj.is_paused = 1;

            obj.current_number = 0;
            
        end
        
        function input(obj, number)
            
            obj.current_number = number;
            
        end
                
        function pause(obj)
           
            obj.stored_time = obj.getElapsed;
            obj.start_time = clock;
            obj.is_paused = 1;
            
        end
        
        function unpause(obj)
           
            obj.start_time = clock;
            obj.is_paused = 0;
            
        end
        
    end
       
    methods % printing
        
        function f = frac(obj)
           
            if isempty(obj.total_number) || obj.total_number<1
                f = 0;
            else
                f = obj.current_number/obj.total_number;
            end
            
        end
        
        function str = show_bar(obj, number)
                        
            if nargin>1 && ~isempty(number)
                obj.current_number = number;
            end
            
            str = sprintf('[%-64s]', repmat('*', 1, ceil(64*obj.frac)));
            
            if nargout<1
                disp(str);
            end
           
        end
        
        function str = show_current_time(obj, number)
                        
            if nargin>1 && ~isempty(number)
                obj.current_number = number;
            end
            
            str = util.text.secs2hms(obj.getElapsed);
            
            if nargout<1
                disp(str);
            end
            
        end   
        
        function str = show_estimate_time(obj, number)
                        
            if nargin>1 && ~isempty(number)
                obj.current_number = number;
            end
            
            str = util.text.secs2hms(obj.getElapsed/obj.frac);
            
            if nargout<1
                disp(str);
            end
            
        end
                        
        function str_out = show(obj, number)
                        
            if nargin>1 && ~isempty(number)
                obj.current_number = number;
            end
            
            str = [obj.show_bar ' ' obj.show_current_time ' / ' obj.show_estimate_time];
            
            if nargout<1 && obj.print_bit
                fprintf('\n');
                disp(str);
                fprintf('\n');
            else
                str_out = str;
            end
            
        end
        
        function c = check(obj, number)
            
            c = number>0 && number~=obj.current_number && mod(number, obj.dividor)==0;
            
        end
        
        function showif(obj, number)
                        
            if isempty(obj.total_number)
                return;
            end
            
            if obj.check(number)
%                 obj.input(number);
                obj.show(number);
            end
            
        end
        
    end   
    
end