classdef WorkingDirectory < handle
% Keeps track of a directory in the filesystem. 
% Lets you change dir, browse, and match files from the directory. 
%
% TEST PROTOCOL: d=util.sys.WorkingDirectory; d.cd('..'); d.match('*.m'); d.pwd

    properties
    
        cwd;
        
    end
    
    properties(Hidden=true)
       
        prev_dir;
        
    end
    
    events
       
        directory_changed;
        
    end
    
    methods % constructor 
       
        function obj = WorkingDirectory(start_at)
            
            if nargin==0 || isempty(start_at)
                obj.cwd = pwd;
            else
                if isa(start_at, 'util.sys.WorkingDirectory')
                    obj.cwd = start_at.cwd;
                elseif ischar(start_at)
                    obj.cwd = start_at;
                else
                    error('input a string or WorkingDirectory object as input to constructor');
                end
            end
            
        end
        
    end
       
    methods % getters
    
        function d = get.cwd(obj)
            
            if strcmp(obj.cwd(end),'/') || strcmp(obj.cwd(end),'/')
                d = obj.cwd;
            elseif strcmp(obj.cwd,'.') || strcmp(obj.cwd, './')
                d = pwd;
            else
                d = [obj.cwd '/'];
            end
            
        end
                  
        function list = ls(obj, directory)
            
            if nargin<2 || isempty(directory)
                directory = obj.cwd;
            end
            
            directory = obj.check_dir(directory);
            
            if nargout==0
                ls(directory);
            else
                list = dir(directory);
                list = {list.name}';
            end
                        
        end
        
        function lh(obj, directory)
            
            if nargin<2 || isempty(directory)
                directory = obj.cwd;
            end
            
            directory = obj.check_dir(directory);
            
            d = dir(directory);
            
            for ii = 1:length(d)
               
                if ~strcmp(d(ii).name,'.') && ~strcmp(d(ii).name,'..')
                    
                    if d(ii).isdir
                        util.text.cprintf('Magenta', '%-60s %-15s \n', d(ii).name, d(ii).date);
                    elseif d(ii).bytes>1024^3
                        fprintf('%-60s %-15s %15f Gbs\n', d(ii).name, d(ii).date, d(ii).bytes/1024^3);
                    elseif d(ii).bytes>1024^2
                        fprintf('%-60s %-15s %15f Mbs\n', d(ii).name, d(ii).date, d(ii).bytes/1024^2);
                    elseif d(ii).bytes>1024
                        fprintf('%-60s %-15s %15f kbs\n', d(ii).name, d(ii).date, d(ii).bytes/1024^1);
                    else
                        fprintf('%-60s %-15s %15f bytess\n', d(ii).name, d(ii).date, d(ii).bytes/1024^0);
                    end
                end
                
            end
                        
        end
        
        function DN_out = dir(obj, directory)
            
            if nargin<2 || isempty(directory)
                directory = obj.cwd;
            end
            
            directory = obj.check_dir(directory);
  
            DN = obj.dir_private(directory);
            
            if nargout==0
                disp(DN);
            else
                DN_out = DN;
            end
            
        end
        
        function DN = dir_private(obj, directory)
                        
            if nargin<2 || isempty(directory)
                directory = obj.cwd;
            end
            
            D = dir(directory);
            N = {D.name}';
            D = {D.isdir}';
            DN = {}; % directory names
            
            for ii = 1:length(N)
                
                if D{ii}==1 && ~strcmp(N{ii}, '.') && ~strcmp(N{ii}, '..')
                    DN{end+1} = N{ii};
                end
                
            end
            
            DN = DN';
            
        end
        
        function list = match(obj, expr, directory)
            
            if nargin<2 
                expr = '';
            end
            
            if nargin<3 || isempty(directory)
                directory = obj.cwd;
            end
                        
            directory = obj.check_dir(directory);
            
            assert(~isempty(directory), ['no such directory: ' directory]);
                        
%             list = util.sys.getfiles(expr, 0, 'glob', directory);

            list = dir(util.text.sa(directory, expr));
            folders = {list.folder}';
            names = {list.name}';
            idx = ~[list.isdir]';
            
            list = strcat(folders(idx), '/', names(idx));
            
        end
        
        function list = files(obj, directory)
                        
            if nargin<2 || isempty(directory)
                directory = obj.cwd;
            end
            
            directory = obj.check_dir(directory);
  
            assert(~isempty(directory), ['no such directory: ' directory]);
            
            FL = obj.files_private(directory);
            
            if nargout==0
                disp(FL);
            else
                list = FL;
            end
            
        end
        
        function list = fullnames(obj, directory)
            
            list = obj.files;
            
            for ii = 1:length(list)
                list{ii} = [obj.cwd list{ii}];
            end
            
        end
        
        function FL = files_private(obj, directory)
                        
            if nargin<2 || isempty(directory)
                directory = obj.cwd;
            end
            
            D = dir(directory);
            N = {D.name}';
            D = {D.isdir}';
            FL = {}; % file list
            
            for ii = 1:length(N)
                
                if D{ii}==0 && ~strcmp(N{ii}, '.') && ~strcmp(N{ii}, '..')
                    FL{end+1} = N{ii};
                end
                
            end
            
            FL = FL';
            
        end
        
        function d = pwd(obj)
            
            d = obj.cwd;
            
            if length(d)>1 && strcmp(d(end), '/')
                d = d(1:end-1);
            end
            
        end
        
        function d = tail(obj)
            
            [~, d] = fileparts(obj.pwd);
            
        end
        
        function d = two_tail(obj)
            
            [a, b] = fileparts(obj.pwd);
            
            [~, c] = fileparts(a);
            
            d = util.text.slash_append(c,b);
            
        end
    end
    
    methods % changing directory
             
        function dir = browse(obj)
           
            dir = uigetdir(obj.pwd);
            
            if ~isempty(dir) && ~isnumeric(dir)
                obj.cd(dir);
            end
            
        end
        
        function CWD = check_dir(obj, next_dir, CWD)
                        
            if nargin<3 || isempty(CWD)
                CWD = obj.cwd;
            end

            if isnumeric(next_dir)
                
                DN = obj.dir_private(CWD);
                
                if next_dir<=length(DN)
                    CWD = [CWD DN{next_dir}];
                    return;
                else
                    CWD = '';
                    return;
%                     error(['current directory ' CWD ' has only ' num2str(length(DN)) ' subdirectories']);
                end
                
            end
            
            if ~strcmp(next_dir(end), '/')
                next_dir = [next_dir '/'];
            end
            
            % if absolute path is used, it is easy...
            if strcmp(next_dir(1),'/') || ~isempty(regexp(next_dir, '^[A-Z]:', 'once'))
                if exist(next_dir,'dir')
                    CWD = next_dir;
                    return;
                else
                    CWD = '';
                    return;
%                     error(['no such directory: ' next_dir]);
                end
            end
            
            % check if "next_dir" is composite
            ind = strfind(next_dir(1:end-1), '/');
            if ~isempty(ind)
                first_dir = next_dir(1:ind);
                second_dir = next_dir(ind+1:end);
            else
                first_dir = next_dir;
                second_dir = '';
            end
            
            if strcmp(first_dir, '..') || strcmp(first_dir, '../') % if the first part is "up"
            
                CWD = [fileparts(CWD(1:end-1)) '/']; % go up one level in CWD
                
                if isempty(CWD)
                    error('cannot go up from here...');
                end
               
            elseif strcmp(first_dir,'.')
                % in this case do nothing...
            else % if the first part is going down the tree
                
                DN = obj.dir_private(CWD);
                
                found = '';
                                
                for ii = 1:length(DN)
                    if strcmp([DN{ii} '/'], first_dir)
                        found = first_dir;
                        break;
                    end
                end % for ii
                
                if isempty(found)
%                     error(['no such directory: ' CWD first_dir]);
                    CWD = '';
                else
                    CWD = [CWD found];
                end
            end
            
            % now we successfully moved to "first_dir"
            if isempty(second_dir)
                return;
            else
                CWD = obj.check_dir(second_dir, CWD);
                return;
            end
            
        end

        function success = cd(obj, move_to)
            
            if isempty(move_to)
                success = 0;
                return;
            end
            
            if isa(move_to, 'util.sys.WorkingDirectory')
                move_to = move_to.pwd;
            end
            
            temp_cwd = obj.check_dir(move_to);
            
            if ~isempty(temp_cwd)
                obj.prev_dir = obj.cwd;
                obj.cwd = temp_cwd;
                
                notify(obj, 'directory_changed');
                
                if nargout>0, success = 1; end
                
            else
                warning(['cannot find directory: ' move_to]);
                
                if nargout>0, success = 0; end
                    
            end
            
        end
        
        function success = smart_cd(obj, str, CWD)
           
            success = 0;
            
            if nargin<3 || isempty(CWD)
                CWD = obj.cwd;
            else
                temp_dir = obj.check_dir(CWD);
                if isempty(CWD)
                    error(['cannot find dir ' CWD]);
                else
                    CWD = temp_dir;
                end
            end
            
            D = obj.dir(CWD);
            
            if exist(util.text.slash_append(CWD, str),'dir')
                obj.cd(util.text.slash_append(CWD, str))
                success = 1;
            else
                
                fnd = [];
                
                for ii = 1:length(D)
                    
                    fnd = strfind(D{ii}, str);
                    
                    if ~isempty(fnd)
                        obj.cd(util.text.slash_append(CWD, D{ii}));
                        success = 1;
                        return;
                    end
                    
                end
                
                if isempty(fnd)
                    disp(['Cannot find expression "' str '" in subdirectories of ' CWD]);
                end
            
            end
            
        end
        
        function here(obj) % change the OBJECT directory to current working dir
           
            obj.cwd = pwd;
            
        end
                
        function there(obj) % change the current working dir to the OBJECT cwd
            
            cd(obj.cwd);
            
        end
    
        function up(obj)
           
            obj.cd('..');
            
        end
        
        function back(obj)
           
            if isempty(obj.prev_dir)
                return;
            end
            
            obj.cd(obj.prev_dir);
            
        end
        
    end
     
    methods % display
        
        function disp(obj)
            
            if isempty(obj)
                builtin('disp', obj);
            else
                display(obj.cwd);
            end
            
        end
                
    end
    
end
    