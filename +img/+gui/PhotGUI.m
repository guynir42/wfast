classdef PhotGUI < handle
    
    properties 
        
        owner@img.Photometry; % link back to containg object

        fig@util.plot.FigHandler;
             
        buttons = {};
        
        font_size = 16;
        edit_font_size = 14;
        small_font_size = 12;
        
        debug_bit = 1;
        
    end
    
    properties % gui stuff
        
        panel_controls;
        panel_aperutres;
        
        panel_info;
        panel_image;
        panel_show;
        axes_image;
        
        panel_close;
        button_close;
        
    end
    
    properties (Hidden=true)
              
        version = 1.00;
        
    end
            
    methods % constructor
       
        function obj = PhotGUI (owner)
            
            % later add other options like copy constructor
            if isa(owner, 'img.Photometry')
                
                if obj.debug_bit, fprintf('PhotGUI  constructor v%4.2f\n', obj.version); end
                
                obj.owner = owner;
                
            else
                error('Input an img.Photometry to constructor of PhotGUI!');
            end
            
        end
        
    end
    
    methods % make/update/check GUI
    
        function make(obj)
            
            import util.plot.GraphicButton;
            import util.plot.GraphicPanel;
            import util.plot.ContrastLimits;
            
            obj.buttons = {};
            
            obj.fig = util.plot.FigHandler('Photometry');
            obj.fig.clear;
            obj.fig.bottom = 5;
            obj.fig.height = 16;
            obj.fig.width = 25;
            movegui(obj.fig.fig, 'center');
            
            %%%%%%%%%%%%%%%%%%% LEFT SIDE %%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            N = 12; % number of buttons on each side
            pos = N;
            
            %%%%%%%%%%% panel controls %%%%%%%%%%%%%%%
            
            % Add buttons using obj.addButton(button_name, var_name='', type='', str1='', str2='', font_size='', split=1, color_on=[], color_off=[])
            
            num_buttons = 6;
            pos = pos-num_buttons;
            obj.panel_controls = GraphicPanel(obj.owner, [0 pos/N 0.2 num_buttons/N], 'controls', 1); % last input is for vertical (default)
            obj.panel_controls.number = num_buttons;
            obj.panel_controls.addButton('button_mex', 'use_mex', 'toggle', 'mex', 'mex', '', [], 'red');
            obj.panel_controls.addButton('button_bg', 'use_backgrounds', 'toggle', 'b/g', 'b/g', '', [], 'red');
            obj.panel_controls.addButton('button_basic', 'use_basic', 'toggle', 'basic', 'basic', '', [], 'red');
            obj.panel_controls.addButton('button_ap', 'use_aperture', 'toggle', 'aperture', 'aperture', '', [], 'red');
            obj.panel_controls.addButton('button_gauss', 'use_gaussian', 'toggle', 'gaussian', 'gaussian', '', [], 'red');
            obj.panel_controls.addButton('input_iter', 'iterations', 'input', 'Niter= ');
            
            obj.panel_controls.make;
            
            %%%%%%%%%%% panel apertures %%%%%%%%%%%%%%
            
            num_buttons = 5;
            pos = pos-num_buttons;
            obj.panel_aperutres = GraphicPanel(obj.owner, [0 pos/N 0.2 num_buttons/N], 'aperture', 1); % last input is for vertical (default)
            obj.panel_aperutres.number = num_buttons;
            obj.panel_aperutres.addButton('button_corners', 'corner_size', 'input', 'corners= ');
            obj.panel_aperutres.addButton('button_aperture', 'aperture', 'input', 'ap= ');
            obj.panel_aperutres.addButton('button_annulus', 'annulus', 'input', 'ann= ');
            obj.panel_aperutres.addButton('button_outer', 'annulus_outer', 'input', 'ann.out= '); 
            obj.panel_aperutres.addButton('button_gauss_sigma', 'gauss_sigma', 'input', 'gauss sig= ');
            
            obj.panel_aperutres.make;
            
            
            %%%%%%%%%%%%%%%%%%% RIGHT SIDE %%%%%%%%%%%%%%%%%%%%%%%%%%
            
            pos = N;
            
            %%%%%%%%%%%%%%%%%%%%%% MIDDLE %%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            obj.panel_info = GraphicPanel(obj.owner, [0.2 (N-1)/N 0.8 1/N], 'info', 0);
            obj.panel_info.number = 6;
            obj.panel_info.addButton('button_flux', 'average_flux', 'info', 'F= ', '', 'small');
            obj.panel_info.addButton('button_bg', 'average_background', 'info', 'B= ', '', 'small');
            obj.panel_info.addButton('button_var', 'average_variance', 'info', 'V= ', '', 'small');
            obj.panel_info.addButton('button_dx', 'average_offset_x', 'info', 'dx= ', '', 'small');
            obj.panel_info.addButton('button_dy', 'average_offset_y', 'info', 'dy= ', '', 'small');
            obj.panel_info.addButton('button_width', 'average_width', 'info', 'W= ', '', 'small');
            obj.panel_info.make;
            
            %%%%%%%%%%% panel image %%%%%%%%%%%%%%%%%%
            
            obj.panel_image = uipanel('Title', '', 'Position', [0.2 1/N 0.8 (N-2)/N]);
            
            %%%%%%%%%%% panel show %%%%%%%%%%%%%%%%%%%
            
            obj.panel_show = GraphicPanel(obj.owner, [0.2 0 0.8 1/N], '', 1);
            obj.panel_show.number = 1;
            obj.panel_show.addButton('button_num_stars', 'show_num_stars', 'input', 'stars= ', '', '', 0.2);
            obj.panel_show.addButton('button_num_frames', 'show_num_frames', 'input', 'frames= ', '', '', 0.2);
            obj.panel_show.addButton('button_picker', '', 'custom', '', '', 'small', 0.6);
            obj.panel_show.make;
            
            %%%%%%%%%%% panel close %%%%%%%%%%%%%%%%%%
            
            obj.panel_close = uipanel('Position', [0 0 0.2 1/N]);
            obj.button_close = GraphicButton(obj.panel_close, [0 0 1 1], obj.owner, '', 'custom', 'CLOSE');
            obj.button_close.Callback = @obj.callback_close;
            
            obj.update;
            
        end
                
        function update(obj,~,~)
                        
            if ~obj.check
                return;
            end
           
            for ii = 1:length(obj.buttons)
                obj.buttons{ii}.update;
            end
            
            obj.owner.plot;
            
        end
                        
        function c = check(obj)
           
            c = ~isempty(obj) && ~isempty(obj.panel_close) && isvalid(obj.panel_close);
            
        end
        
    end
                
    methods % callbacks
        
        function callback_close(obj, ~, ~)
           
            if obj.debug_bit, disp('callback: close'); end
            
            delete(obj.fig.fig);
            
        end
        
    end
    
end