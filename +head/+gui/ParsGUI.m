classdef ParsGUI < handle
    
    properties 
        
        owner@head.Parameters; % link back to containg object

        fig@util.plot.FigHandler;
             
        buttons = {};
        buttons_stars = {};
        
        font_size = 16;
        edit_font_size = 14;
        small_font_size = 12;
        
        debug_bit = 1;
        
    end
    
    properties % gui stuff
        
        panel_target;
        panel_instrument;
        panel_stars;
        vector_star_buttons;
        button_add_star;
        button_clear_stars;
            
        panel_close;
        button_close;
            
    end
    
    properties (Hidden=true)
              
        version = 1.00;
        
    end
            
    methods % constructor
       
        function obj = ParsGUI(owner)
            
            % later add other options like copy constructor
            if isa(owner, 'head.Parameters')
                
                if obj.debug_bit, fprintf('ParsGUI constructor v%4.2f\n', obj.version); end
                
                obj.owner = owner;
                                
            else
                error('Input a head.Parameters object to constructor of ParsGUI!');
            end
            
        end
        
    end
    
    methods % make/update/check GUI
    
        function make(obj)
            
            import util.plot.GraphicButton;
            import util.plot.GraphicPanel;
            import util.plot.ContrastLimits;
            
            obj.buttons = {};
            obj.buttons_stars = {};
            
            if isempty(obj.fig)
                obj.fig = util.plot.FigHandler('Parameters');
            end
            
            obj.fig.reset;
            obj.fig.left = 5;
            obj.fig.bottom = 5;
            obj.fig.height = 16;
            obj.fig.width = 10;
%             obj.fig.name = '...';
            movegui(obj.fig.fig, 'center');
            
            %%%%%%%%%%% panel target %%%%%%%%%%%%%%%%%%
            
            obj.panel_target = GraphicPanel(obj.owner, [0, 0.7, 0.5, 0.3], 'Target');
            obj.panel_target.addButton('button_reset', 'resetTarget', 'push', 'RESET');
            obj.panel_target.addButton('button_name', 'target_name', 'input',' ', '', '', 0.8);
            obj.panel_target.addButton('button_push_name', 'name_plus_plus', 'push', '++','','', 0.2);
            obj.panel_target.addButton('button_RA', 'RA', 'input', 'RA: ');
            obj.panel_target.addButton('button_DE', 'DE', 'input', 'DE: ');
            obj.panel_target.number = 5;
            obj.panel_target.make;
            obj.panel_target.button_name.Callback = @obj.callback_target_name;
            
            %%%%%%%%%%% panel instrument %%%%%%%%%%%%%%
            
            obj.panel_instrument = GraphicPanel(obj.owner, [0.5 0.7 0.5 0.3], 'Instrument');
            obj.panel_instrument.addButton('button_filter', 'filter_name', 'input', ' ');
            obj.panel_instrument.addButton('button_aperture', 'aperture', 'input', 'D= ', ' cm');
            obj.panel_instrument.addButton('button_f_number', 'f_number', 'input', 'f/#= ');
            obj.panel_instrument.number = 5;
            obj.panel_instrument.make;
            obj.panel_instrument.button_filter.Callback = @obj.callback_filter_name;
            
            %%%%%%%%%%% panel stars %%%%%%%%%%%%%%%%%%%
            
            obj.panel_stars = uipanel('Title', 'Stars', 'Position', [0 0.1 1 0.6]);
            
            obj.fillStarsPanel;
            
            %%%%%%%%%%% panel close %%%%%%%%%%%%%%%%%%%
            
            obj.panel_close = uipanel('Position', [0 0 1 0.1]);
            
            obj.button_close = GraphicButton(obj.panel_close, [0 0 1 1], obj.owner, '', 'custom', 'CLOSE');
            obj.button_close.Callback = @obj.callback_close;
            
            obj.update;
            
        end
        
        function fillStarsPanel(obj)
           
            import util.plot.GraphicButton;
            
            N = 8;
            
            % get rid of previous buttons...
            for ii = 1:length(obj.vector_star_buttons)
                delete(obj.vector_star_buttons(ii).control);
            end
            
            if ~isempty(obj.button_add_star), delete(obj.button_add_star.control); end
            if ~isempty(obj.button_add_star), delete(obj.button_clear_stars.control); end
            
            obj.vector_star_buttons = GraphicButton.empty;
            obj.buttons_stars = {};
            
            % make new buttons
            for ii = 1:length(obj.owner.stars)
                
                str = obj.owner.stars(ii).printout;
                obj.vector_star_buttons(ii) = GraphicButton(obj.panel_stars, [0 (N-ii)/N 1 1/N], obj.owner, '', 'custom', str, '','', '','buttons_stars');
                obj.vector_star_buttons(ii).control.UserData = ii;
                obj.vector_star_buttons(ii).Callback = @obj.callback_makeStarGUI;
                if ii>N-2, break; end
                
            end
            
            Nstars = length(obj.vector_star_buttons);
            
            obj.button_add_star = GraphicButton(obj.panel_stars, [0 (8-Nstars-1)/8 1 1/N], obj.owner, '', 'custom', 'ADD STAR', '','','','buttons_stars');
            obj.button_add_star.Callback = @obj.callback_add_star;
            obj.button_clear_stars = GraphicButton(obj.panel_stars, [0 0 1 1/N], obj.owner, '', 'custom', 'CLEAR STARS', '','','','buttons_stars');
            obj.button_clear_stars.Callback = @obj.callback_clear_stars;
            
        end
        
        function update(obj,~,~)
                        
            if ~obj.check
                return;
            end
           
            for ii = 1:length(obj.buttons)
                obj.buttons{ii}.update;
            end
            
            for ii = 1:length(obj.buttons_stars)-2
                
                obj.buttons_stars{ii}.String = obj.owner.stars(ii).printout;
                obj.buttons_stars{ii}.update;
                
            end
            
        end
                       
        function c = check(obj)
           
            c = ~isempty(obj) && ~isempty(obj.panel_target) && isvalid(obj.panel_target);
            
        end
        
    end
                
    methods % callbacks
        
        function callback_target_name(obj, hndl, ~)
            
            if obj.debug_bit, disp(['hndl.String= ' hndl.String]); end
                        
            obj.owner.target_name = strip(hndl.String);
            
            obj.update;
            
        end
        
        function callback_filter_name(obj, hndl, ~)
            
            if obj.debug_bit, disp(['hndl.String= ' hndl.String]); end
                        
            obj.owner.filter_name = strip(hndl.String);
            
            obj.update;
            
        end
        
        function callback_makeStarGUI(obj, hndl, ~)
           
            value = hndl.UserData;
            if obj.debug_bit, disp(['callback: makeStarGUI... value= ' num2str(value)]); end
            
            obj.owner.stars(value).makeGUI;
            obj.owner.stars(value).gui.owner = obj.owner.stars(value); % make the GUI open up on the right star
            obj.owner.stars(value).gui.update;
            
            obj.update;
            
        end
        
        function callback_add_star(obj, ~, ~)
            
            obj.owner.addStar;
            obj.owner.stars(end).makeGUI;
            obj.owner.stars(end).gui.owner = obj.owner.stars(end);
            obj.owner.stars(end).gui.update;
            obj.fillStarsPanel;
            
            obj.update;
            
        end
        
        function callback_clear_stars(obj, ~, ~)
            
            obj.owner.stars = head.Star.empty;
            obj.fillStarsPanel;
            obj.update;
            
        end
        
        function callback_close(obj, ~, ~)
           
            if obj.debug_bit, disp('callback: close'); end
            
            delete(obj.fig.fig);
            
        end
        
    end
    
end