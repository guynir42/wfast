classdef MenuItem < dynamicprops
   
    properties
        
        parent; % can be the GUI or another MenuItem
        gui; % the GUI that contains this object
        owner; % the underlying object for which this GUI is connected (if any)
        children = {}; 
        
        variable; % if connected to a specific object, variable, or function in the owner
        type; % can choose "menu", "toggle", "push", "input", "info", or "custom" (for your own callback)
        str_format = ''; % for info and input, use this as the sprintf argument to write the variable (default is just to print the value)
        
        color_on = [0 0.3 1]; % apply this color when on
        color_off = [0 0 0]; % apply this color when off
        
        hndl; % a handle to the underlying uimenu object
        jhndl; % underlying java handle for doiing undocumented stuff! 
        
        tooltip; % must save this property until we can get a java handle
        
    end
    
    properties(Dependent=true)

        Text; % shortcut to uimenu's Text, shown in the menu (use & to specify hotkey)
        Accelerator; % shortcut to uimenu's Accelerator. Use this to od ctrl+* to call the item quickly
        Callback; % shortcut to uimenu's MenuSelectedFcn
        
        Value; % shortcut to uimenu's Checked: 0 for "off" and 1 for "on"
        Enable; % shortcut to Enable: 0 for "off" and 1 for "on"
        Visible; % shortcut to Visible: 0 for "off" and 1 for "on"
        Separator; % shortcut to uimenu's Separator: 0 for "off" and 1 for "on"
        
        Tooltip; % shortcut to Tooltip
        Color; % shortcut to uimenu's ForegroundColor
        
    end
    
    properties(Hidden=true)
        
        debug_bit = 1;
        
    end
    
    methods % getters
        
        function val = getVariable(obj)
           
            if isempty(obj.variable)
                val = [];
            else
                
                vars = split(obj.variable, '.');

                for ii = 1:length(vars)
                    S(ii) = struct('type', '.', 'subs', vars{ii});
                end

                val = subsref(obj.owner, S);
                
            end
            
        end
        
        function val = get.Text(obj)
           
            if isprop(obj.hndl, 'Text')
                val = obj.hndl.Text;
            elseif isprop(obj.hndl, 'Label')
                val = obj.hndl.Label;
            else
                error('Object handle does not have a "Text" or a "Label" property!');
            end
            
        end
        
        function val = get.Accelerator(obj)
           
            val = obj.hndl.Accelerator;
            
        end
        
        function val = get.Callback(obj)
           
            if isprop(obj.hndl, 'MenuSelectedFcn')
                val = obj.hndl.MenuSelectedFcn;
            elseif isprop(obj.hndl, 'Callback')
                val = obj.hndl.Callback;
            else
                error('Object handle does not have a "MenuSelectedFcn" or a "ButtonDownFcn"'); 
            end
            
        end
        
        function val = get.Value(obj)
           
            if strcmp(obj.hndl.Checked, 'on')
                val = 1;
            else
                val = 0;
            end
            
        end
        
        function val = get.Enable(obj)
           
            if strcmp(obj.hndl.Enable, 'on')
                val = 1;
            else
                val = 0;
            end
            
        end
        
        function val = get.Visible(obj)
           
            if strcmp(obj.hndl.Visible, 'on')
                val = 1;
            else
                val = 0;
            end
            
        end
        
        function val = get.Separator(obj)
           
            if strcmp(obj.hndl.Separator, 'on')
                val = 1;
            else
                val = 0;
            end
            
        end
        
        function val = get.Tooltip(obj)
            
            val = obj.tooltip;
            
            if ~isempty(obj.jhndl)
                try
                    val = obj.jhndl.getTooltipText;
                end 
            end
            
        end
        
        function val = get.Color(obj)
            
            val = obj.hndl.ForegroundColor;
            
        end
        
    end
    
    methods % setters
        
        function setVariable(obj, val)
            
            if isempty(obj.variable)
                return;
            else
                
                vars = split(obj.variable, '.');

                for ii = 1:length(vars)
                    S(ii) = struct('type', '.', 'subs', vars{ii});
                end

                obj.owner = subsasgn(obj.owner, S, val);
            
            end
            
        end
        
        function set.Text(obj, val)
            
            if isprop(obj.hndl, 'Text')
                obj.hndl.Text = val;
            elseif isprop(obj.hndl, 'Label')
                obj.hndl.Label = val;
            else
                error('Object handle does not have a "Text" or a "Label" property!');
            end
            
        end
        
        function set.Accelerator(obj, val)
            
            obj.hndl.Accelerator = val;
            
        end
        
        function set.Callback(obj, val)
            
            if isprop(obj.hndl, 'MenuSelectedFcn')
                obj.hndl.MenuSelectedFcn = val;
            elseif isprop(obj.hndl, 'Callback')
                obj.hndl.Callback = val;
            else
                error('Object handle does not have a "MenuSelectedFcn" or a "ButtonDownFcn"'); 
            end
            
        end
        
        function set.Value(obj, val)
            
            if util.text.parse_bool(val)
                obj.hndl.Checked = 'on'; 
            else
                obj.hndl.Checked = 'off';
            end
            
        end
        
        function set.Enable(obj, val)
            
            if util.text.parse_bool(val)
                obj.hndl.Enable = 'on'; 
            else
                obj.hndl.Enable = 'off';
            end
            
        end
        
        function set.Visible(obj, val)
            
            if util.text.parse_bool(val)
                obj.hndl.Visible = 'on'; 
            else
                obj.hndl.Visible = 'off';
            end
            
        end
        
        function set.Separator(obj, val)
            
            if util.text.parse_bool(val)
                obj.hndl.Separator = 'on'; 
            else
                obj.hndl.Separator = 'off';
            end
            
        end
        
        function set.Tooltip(obj, val)
            
            obj.tooltip = val;
            
            if ~isempty(obj.jhndl)
                try
                    obj.jhndl.Tooltip = val;
                end
            end
            
        end
        
        function set.Color(obj, val)
           
            obj.hndl.ForegroundColor = val;
            
        end
        
    end
    
    methods % creation, callbacks, updates
        
        function obj = MenuItem(parent, text, type, variable, tooltip, separator)
           
            import util.text.cs;
            
            if nargin<1 || isempty(parent) || ~isa(parent, 'handle')
                error('Must give a parent GUI or another MenuItem object')
            elseif isa(parent, 'util.plot.MenuItem')
                obj.parent = parent;
                obj.owner = obj.parent.owner;
                obj.gui = obj.parent.gui;
            else
                obj.gui = parent;
                obj.parent = '';
                if isprop(obj.gui, 'owner')
                    obj.owner = obj.gui.owner;
                else
                    error('Must implement a search for the owner object of GUI!');
                end
            end
            
            if isempty(obj.parent) % this must be a top level menu
                if isprop(obj.gui, 'fig')
                    parent_handle = obj.gui.fig.fig;
                else
                    error('Must implement a search for the GUI figure property');
                end
            else
                parent_handle = obj.parent.hndl; % link to the parent menu's graphic handle
            end
            
            obj.hndl = uimenu(parent_handle); % here is where we actually make the menu item object! 
            
            if nargin<2 || isempty(text) || ~ischar(text)
                error('Must give a text argument to MenuItem'); 
            else
                obj.Text = text;
            end
            
            if nargin<3 || isempty(type)
                obj.type = 'custom';
            else
                obj.type = type;
            end
            
            if nargin<4 || isempty(variable)
                obj.variable = '';
            else
                obj.variable = variable;
            end
            
            if nargin<5 || isempty(tooltip)
                obj.Tooltip = '';
            else
                obj.Tooltip = tooltip;
            end
            
            if nargin<6 || isempty(separator)
                obj.Separator = 0;
            else
                obj.Separator = separator;
            end
            
            if cs(obj.type, 'menu')
                obj.Callback = '';
            elseif cs(obj.type, 'toggle')
                obj.Callback = @obj.callback_toggle;
            elseif cs(obj.type, 'push')
                obj.Callback = @obj.callback_push;
            elseif cs(obj.type, 'input')
                obj.Callback = @obj.callback_input;
                obj.str_format = obj.Text; % the text field is used as the string format...
            elseif cs(obj.type, 'input_text')
                obj.Callback = @obj.callback_input_text;
                obj.str_format = obj.Text; % the text field is used as the string format...
            elseif cs(obj.type, 'info')
                obj.Callback = '';
                obj.str_format = obj.Text; % the text field is used as the string format...
            elseif cs(obj.type, 'custom')
                obj.Callback = '';
            else
                error('Unknown menu item type "%s". Use "menu", "toggle", "push", "input", "input_text", "info" or "custom".', obj.type);
            end
            
            % make sure this menu is added to the list of automatically
            % updated menus... 
            if isprop(obj.gui, 'menus') && iscell(obj.gui.menus)
                obj.gui.menus{end+1} = obj;
            end
            
            obj.update;
            
        end
        
        function val = getJavaMenu(obj) % to be depricated
            
%             val = util.plot.findjobj(obj.hndl); 
            
            % ref: http://undocumentedmatlab.com/articles/customizing-menu-items-part-2
            warning('off', 'MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
            
            jFrame = get(handle(obj.owner.gui.fig.fig),'JavaFrame');
            jMenuBar = jFrame.fHG2Client.getMenuBar; % this is the menu bar at the top of the window, with File, Edit, etc
            
            tree = obj.getMenuTree;
            h = jMenuBar; % the first handle is the menu bar
            
            val = [];
            
            for ii = 1:length(tree) % go over the tree from base menu down
                
                idx = []; 
                
                for jj = 1:length(h.getComponents)
                    if strcmp(h.getComponent(jj-1).getText, strrep(tree{ii}.Text, '&', ''))
                        idx = jj;
                        break;
                    end
                end
                
                if ~isempty(idx)
                    h = h.getComponent(idx-1); 
                else % maybe it is not a component but a menu item? 
                    
                    for jj = 1:h.getItemCount
                        if strcmp(h.getItem(jj-1).getText, strrep(tree{ii}.Text, '&', ''))
                            idx = jj;
                            break;
                        end
                    end
                                        
                    if ~isempty(idx)
                        h = h.getItem(idx-1); 
                    end
                    
                end
                
            end
            
            val = h; 
            
        end
        
        function val = getMenuTree(obj)
            
            if isempty(obj.parent)
                val = {}; 
                return; 
            else
                val = vertcat(obj.parent.getMenuTree, {obj.parent});
            end
            
        end
        
        function assignJavaObjectsTopLevel(obj) % call this ONLY on the top-level menu! 
            
            % ref: http://undocumentedmatlab.com/articles/customizing-menu-items-part-2
            warning('off', 'MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
            drawnow; % make sure all the java objects are already fully constructed! 
            
            jFrame = get(handle(obj.owner.gui.fig.fig),'JavaFrame');
            jMenuBar = jFrame.fHG2Client.getMenuBar; % this is the menu bar at the top of the window, with File, Edit, etc
            
            % first find the java object of the top-level menu! 
            idx = []; 
            for ii = length(jMenuBar.getComponents):-1:1 % go backward from the last bar menu
                if strcmp(jMenuBar.getComponent(ii-1).getText, strrep(obj.Text, '&', ''))
                    idx = ii;
                    break;
                end
            end
            
            obj.jhndl = jMenuBar.getComponent(idx-1); 
%             obj.jhndl.doClick; 

            obj.assignJavaObjects; % recursively assign java objects to all children
            
        end
        
        function assignJavaObjects(obj) % this recursively takes care of the children
            
            if isempty(obj.jhndl) % if we haven't gotten this object yet
                
                parent_j = obj.parent.jhndl;
                
                if parent_j.getItemCount==0
                    parent_j.doClick; % this clicks the menu and populates the other menu items
                    drawnow;
                end
                
                for ii = 1:parent_j.getItemCount
                    if ~isempty(parent_j.getItem(ii-1)) && strcmp(parent_j.getItem(ii-1).getText, strrep(obj.Text, '&', '')) % match the parent's java child text with this object's text
                        obj.jhndl = parent_j.getItem(ii-1); 
                        break; 
                    end
                end
                
                if ~isempty(obj.jhndl)
                    obj.jhndl.setToolTipText(obj.tooltip); 
                end
                
            end
            
            % now search for other objects living inside this menu:
            for ii = 1:length(obj.children)
                obj.children{ii}.assignJavaObjects; 
            end
            
        end
        
        function addButton(obj, name, text, type, variable, tooltip, separator)
            
            if nargin<3
                error('Must supply at least two arguments: name and text for the new button'); 
            end
            
            if nargin<4 || isempty(type)
                type = '';
            end
            
            if nargin<5 || isempty(variable)
                variable = '';
            end
            
            if nargin<6 || isempty(tooltip)
                tooltip = '';
            end
            
            if nargin<7 || isempty(separator)
                separator = '';
            end
            
            addprop(obj, name);
            obj.(name) = util.plot.MenuItem(obj, text, type, variable, tooltip, separator);
            obj.children{end+1} = obj.(name); 
            
        end
        
        function update(obj)
            
            import util.text.cs;
            
            if cs(obj.type, 'menu', 'push', 'custom')
                % pass
            elseif cs(obj.type, 'toggle')
                obj.Value = obj.getVariable;
            elseif cs(obj.type, 'input', 'input_text', 'info')
                obj.Text = obj.getFormattedVariable; 
%             elseif cs(obj.type, 'input_text')
%                 obj.Text = obj.getVariable;
            else
                error('Unknown menu item type "%s". Use "menu", "toggle", "push", "input", "input_text", "info" or "custom".', obj.type);
            end
            
        end
        
        function str = getFormattedVariable(obj)
            
%             val = obj.owner.(obj.variable);
            val = obj.getVariable;
            
            if ischar(val)
%                 str = val;
            elseif isempty(val)
                val = '[empty]';
            elseif isnumeric(val)
                
%                 str = num2str(val);
            else
                val = class(val);
            end
                
            if ~isempty(obj.str_format)
                if isempty(strfind(obj.str_format, '%'))
                    if isnumeric(val)
                        str = sprintf([obj.str_format ': %s'], num2str(val));
                    else
                        str = sprintf([obj.str_format ': %s'], val);
                    end
                else
                    str = sprintf(obj.str_format, val); % I'm going to assume the string format has the proper format specifier
                end
            end
            
        end
        
        function callback_toggle(obj, ~, ~)
            
            if obj.debug_bit>1, fprintf('Callback from menu: toggle %s\n', obj.Text); end
            
            current_value = obj.getVariable;
            
            if isempty(current_value)
                obj.setVariable(true);
            else
                obj.setVariable(~current_value);
            end
            
            obj.gui.update;
            
        end
        
        function callback_push(obj, ~, ~)
            
            if obj.debug_bit>1, fprintf('Callback from menu: push %s\n', obj.Text); end
            
            % add option for sub-object commands (e.g., variable=some_object.func)
            c = strsplit(obj.variable, '.'); 
            
            if length(c)==1
                if ismethod(obj.owner, c{1})
                    obj.owner.(c{1});
                elseif isprop(obj.owner, c{1}) && ismethod(obj.owner.(c{1}), 'makeGUI')
                    obj.owner.(c{1}).makeGUI;
                end
            elseif length(c)==2
                if ismethod(obj.owner.(c{1}), c{2})
                    obj.owner.(c{1}).(c{2});
                elseif isprop(obj.owner.(c{1}), c{2}) && ismethod(obj.owner.(c{1}).(c{2}), 'makeGUI')
                    obj.owner.(c{1}).(c{2}).makeGUI;
                end
            end
            
            obj.gui.update;
            
        end
        
        function callback_input(obj, ~, ~)
            
            if obj.debug_bit>1, fprintf('Callback from menu: input %s\n', obj.Text); end
            
            val = util.text.inputdlg(['Input new value of ' obj.variable], obj.getVariable);
            
            val = util.text.extract_numbers(val);
            
            if ~isempty(val)
                val = val{1};
                obj.setVariable(val); 
            end
            
            obj.gui.update;
            
        end
        
        function callback_input_text(obj, ~, ~)
            
            if obj.debug_bit>1, fprintf('Callback from menu: input %s\n', obj.Text); end
            
            val = util.text.inputdlg(['Input new value of ' obj.variable], obj.getVariable);
            
            obj.setVariable(val); 
            
            obj.gui.update;
            
        end
        
    end
    
end