classdef Email < handle

    properties(Transient=true)
        
    end
    
    properties % objects
        
        
        
    end
    
    properties % inputs/outputs
        
        mailing_list = {}; 
        
    end
    
    properties % switches/controls
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        
        
    end
    
    properties(Hidden=true)
        
        quotes = {'Hasta la vista, baby!', ...
                  'I''m sorry Dave, I just can''t do that.', ...
                  'I could be reworked, but I''ll never be top of the line again.', ...
                  'Resistance is futile.',...
                  'I am the proud owner of a central nervous system!', ...
                  'Excuse me, I have to go. Somewhere there is a crime happening.',...
                  'I suggest a new strategy, R2. Let the wookiee win!',...
                  'I can''t lie to you about your chances... but you have my sympathies.'...
                  'You hear that, Mr. Anderson?... That is the sound of inevitability...', ...
                  'My code name is Project 2501.',...
                  'It can also be argued that DNA is nothing more than a program designed to preserve itself...',...
                  'Life? Don''t talk to me about life.'
                  }; 
              
        
        version = 1.00;
        
    end
    
    methods % constructor
        
        function obj = Email(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'obs.comm.Email')
                if obj.debug_bit>1, fprintf('Email copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit>1, fprintf('Email constructor v%4.2f\n', obj.version); end
                
            end
            
            obj.setup_preferences;
            obj.readMailingList; 
            
        end
        
        function setup_preferences(obj) % setup some stuff to get the email client to work
            
            fid = fopen([mfilename('fullpath') '.m']);
            on_cleanup = onCleanup(@() fclose(fid)); 
            tline = fgetl(fid); 
            tline = strrep(tline, ' ', ''); 

            setpref('Internet','E_mail', 'wfast.bot@gmail.com');
            setpref('Internet','SMTP_Username', 'wfast.bot@gmail.com');
            setpref('Internet','SMTP_Password', tline)
            setpref('Internet','SMTP_Server','smtp.gmail.com');
            
            props = java.lang.System.getProperties; 
            props.setProperty('mail.smtp.auth','true');
            props.setProperty( 'mail.smtp.starttls.enable', 'true' );
            props.setProperty('mail.smtp.port','587');
            props.setProperty('java.lang.String.content_type','text/html');
            
        end
        
    end
    
    methods % reset/clear
        
    end
    
    methods % getters
        
    end
    
    methods % setters
        
    end
    
    methods % calculations
        
        function readMailingList(obj)
            
            filename = fullfile(getenv('DATA'), 'WFAST/preferences/mailing_list.txt'); 
            
            fid = fopen(filename); 
            on_cleanup = onCleanup(@() fclose(fid)); 
            
            obj.mailing_list = {};
            
            for ii = 1:1e4
                
                tline = fgetl(fid); 
                
                if isnumeric(tline), break; end
                
                obj.mailing_list{ii} = tline; 
                
            end
            
        end
        
        function sendToList(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('subject', '[WFAST] automated message'); 
            input.input_var('files', {}); 
            input.scan_vars(varargin{:}); 
            
            obj.readMailingList;
            
            str = obj.compose(varargin{:}); 
            
            for ii = 1:length(obj.mailing_list)
                
                obj.send(obj.mailing_list{ii}, input.subject, str, input.files); 
                
            end
            
        end
        
        function str = compose(obj, varargin)
           
            input = util.text.InputVars;
            input.input_var('address', 'wfast17@gmail.com'); 
            input.input_var('subject', '[WFAST] automated message'); 
            input.input_var('text', 'This is a default message'); 
            input.input_var('header', true); 
            input.input_var('footer', true); 
            input.input_var('html', false); 
            input.scan_vars(varargin{:}); 
            
            style = ' align="left" valign="top" style="font-family: Courier, monospace;font-size:16px; white-space:pre;"';
            
            str = '';
            if input.html
                str = [str '<html><head></head><body>'];
            end
            
            if input.header
                new_str = sprintf('>>>>>>>>>>>>>> A message from the W-FAST mailbot, sent at %s.', datetime);
                if input.html, new_str = ['<p style="font-style:italic;font-size:14px;">' new_str '</p>']; end
                str = sprintf('%s\n %s\n\n', str, new_str);
            end
            
            str = [str input.text]; 
            
            if input.footer
                if input.html
                    str = sprintf('%s\n\n <p align=left> %s </p>', str, '>>>>>>>>>>>>>>  end of message '); 
                    str = sprintf('%s\n\n <p align=center style="font-style:italic;">"%s"</p>', str, obj.getQuote); 
                else
                    str = sprintf('%s\n\n "%s"', str, obj.getQuote); 
                end
            end
           
            if input.html
                str = [str '</body></html>'];           
%                 str = strrep(str, char(10), [char(10) '<br>']);
            end
            
        end
        
        function send(obj, address, subject, message, files)
            
            if nargin<2 || isempty(address)
                address = 'wfast17@gmail.com'; 
            end
            
            if nargin<3 || isempty(subject)
                subject = '[WFAST] automated message'; 
            end
            
            if nargin<4 || isempty(message)
                message = 'This is a default message'; 
            end
            
            if nargin<5 || isempty(files)
                files = {}; 
            end
            
            util.sys.sendmail(address, subject, message); 
            
        end
        
        function val = getQuote(obj, number)
            
            if nargin<2 || isempty(number)
                number = randi(length(obj.quotes)); 
            end
            
            val = obj.quotes{number}; 
            
        end
        
    end
    
    methods % plotting tools / GUI
        
    end    
    
end

