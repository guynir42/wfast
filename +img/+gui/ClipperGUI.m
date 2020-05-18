classdef ClipperGUI < handle
   
    properties % objects and other controls
    
        clip@img.Clipper;
        fig@util.plot.FigHandler;
        
        buttons = {};
        
        font_size = 18;
        edit_font_size = 14;
        small_font_size = 12;
        
        debug_bit = 1;
        
    end
    
    properties % gui objects
        
        panel_controls;
        input_cut_size
        button_mex_cut;
        
        panel_adjust;
        button_use_lock_adjust;
        button_use_moments;
                
        panel_find_stars;
        input_num_stars;
        input_filter_size;
        input_filter_sigma;
        input_avoid_edges;
        
        panel_display;
        input_number_cuts;
        
        panel_close;
        button_close;
        
        panel_subframe;
        button_choose_subframe;
        input_subframe_size;
        
        panel_positions;
        button_positions;
        
        panel_plot;
        axes_plot;
        
    end
    
    properties(Hidden=true)
       
        version = 1.01;
        
    end
    
    methods % constructor
        
        function obj = ClipperGUI(clip_obj)
           
            assert(isa(clip_obj, 'img.Clipper'), 'must input a Clipper object to construct a ClipperGUI');
            
            if obj.debug_bit>1, fprintf('ClipperGUI constructor v%4.2f\n', obj.version); end
            
            obj.clip = clip_obj;
                        
        end
        
    end
    
    methods % make and update GUI
        
        function makeGUI(obj)
                      
            import util.plot.GraphicButton;
            
            if isempty(obj.fig)
                obj.fig = util.plot.FigHandler('clipper GUI');
            end
            
            obj.buttons = {};
            
            obj.fig.reset;
            obj.fig.name = 'Clipper';
            
            %%%%%%%%%%%%%%%%%% panel controls %%%%%%%%%%%%%%%%%%%%%%
            
            obj.panel_controls = uipanel('Title', 'controls', 'Units', 'Normalized', 'Position', [0 0.8 0.2 0.2]);
            
            obj.input_cut_size = uicontrol(obj.panel_controls, 'Style', 'edit',...
                'Units', 'Normalized', 'Position', [0 0/2 1 1/2],...
                'Callback', @obj.callback_cut_size);
            
            obj.button_mex_cut = uicontrol(obj.panel_controls, 'Style', 'pushbutton',...
                'Units', 'Normalized', 'Position', [0 1/2 1 1/2],...
                'Callback', @obj.callback_mex_cut);
            
            %%%%%%%%%%%%%%%%%% panel adjust %%%%%%%%%%%%%%%%%%%%%%%%
            
            obj.panel_adjust = uipanel('Title', 'adjust', 'Units', 'Normalized', 'Position', [0 0.6 0.2 0.2]);
            
            obj.button_use_lock_adjust = uicontrol(obj.panel_adjust, 'Style', 'pushbutton',...
                'Units', 'Normalized', 'Position', [0 0/2 1 1/2],...
                'Callback', @obj.callback_use_lock_adjust);
            
            obj.button_use_moments = uicontrol(obj.panel_adjust, 'Style', 'pushbutton',...
                'Units', 'Normalized', 'Position', [0 1/2 1 1/2],...
                'Callback', @obj.callback_use_moments);
            
            %%%%%%%%%%%%%%%%%% panel find stars %%%%%%%%%%%%%%%%%%%%
            
            obj.panel_find_stars = uipanel('Title', 'find stars', 'Units', 'Normalized', 'Position', [0 0.2 0.2 0.4]);
            
            obj.input_num_stars = uicontrol(obj.panel_find_stars, 'Style', 'edit',...
                'Units', 'Normalized', 'Position', [0 0/4 1 1/4],...
                'Callback', @obj.callback_num_stars);

            obj.input_filter_size = uicontrol(obj.panel_find_stars, 'Style', 'edit',...
                'Units', 'Normalized', 'Position', [0 1/4 1 1/4],...
                'Callback', @obj.callback_filter_size);

            obj.input_filter_sigma = uicontrol(obj.panel_find_stars, 'Style', 'edit',...
                'Units', 'Normalized', 'Position', [0 2/4 1 1/4],...
                'Callback', @obj.callback_filter_sigma);

            obj.input_avoid_edges = uicontrol(obj.panel_find_stars, 'Style', 'edit',...
                'Units', 'Normalized', 'Position', [0 3/4 1 1/4],...
                'Callback', @obj.callback_avoid_edges);
                        
            %%%%%%%%%%%%%%%%%% panel display %%%%%%%%%%%%%%%%%%%%%%%
            
            obj.panel_display = uipanel('Title', 'display', 'Units', 'Normalized', 'Position', [0.8 0.8 0.2 0.2]);
            
            num = 2;
            
            obj.input_number_cuts = GraphicButton(obj.panel_display, [0 1/num 1 1/num], obj.clip, 'number_cuts_display', 'input', 'Ndisp= ');
            
            %%%%%%%%%%%%%%%%%% panel subframe %%%%%%%%%%%%%%%%%%%%%%
            
            obj.panel_subframe = uipanel('Title', 'subframe', 'Units', 'Normalized', 'Position', [0.8 0.6 0.2 0.2]);
            
            num = 2;
            
            obj.button_choose_subframe = GraphicButton(obj.panel_subframe, [0 1/num 1 1/num], obj.clip, 'chooseSubframe','push','choose');
            obj.input_subframe_size = GraphicButton(obj.panel_subframe, [0 0/num 1 1/num], obj.clip, 'subframe_size', 'input', 'size= ');
            
            %%%%%%%%%%%%%%%%%% panel positions %%%%%%%%%%%%%%%%%%%%%
            
            obj.panel_positions = uipanel('Title', 'positions', 'Units', 'Normalized', 'Position', [0.8 0 0.2 0.6]);
            
            obj.button_positions = uicontrol(obj.panel_positions, 'Style', 'pushbutton',...
                'Units', 'Normalized', 'Position', [0 0 1 1],...
                'Callback', @obj.update);
            
            %%%%%%%%%%%%%%%%%% panel close %%%%%%%%%%%%%%%%%%%%%%%%%
            
            obj.panel_close = uipanel('Title', '', 'Units', 'Normalized', 'Position', [0 0.0 0.2 0.2]);
            
            obj.button_close = uicontrol(obj.panel_close, 'Style', 'pushbutton',...
                'Units', 'Normalized', 'Position', [0 0 1 1],...
                'Callback', @obj.callback_close, 'String', 'CLOSE');
            
            %%%%%%%%%%%%%%%%%% panel plot %%%%%%%%%%%%%%%%%%%%%%%%%%
            
            obj.panel_plot = uipanel('Title', '', 'Units', 'Normalized', 'Position', [0.2 0 0.6 1]);
  
            obj.makeAxes;
            
%             title('drifts in cut positions');
%             xlabel('frame number');
%             ylabel('pixel drift');

            obj.update;
            obj.clip.show;
            
        end
        
        function makeAxes(obj)
            
                      
            delete(obj.panel_plot.Children);
            obj.axes_plot = {};
            Nrows = ceil(sqrt(obj.clip.number_cuts_display));
            Ncols = Nrows;
            
            for ii = 1:obj.clip.number_cuts_display
                
                x = mod(ii-1, Nrows);
                y = floor((ii-1)/Nrows);
                
                obj.axes_plot{ii} = axes('Parent', obj.panel_plot, 'Position', [x/Ncols y/Nrows 1/Ncols 1/Nrows]);
                obj.axes_plot{ii}.XTick = [];
                obj.axes_plot{ii}.YTick = [];
%                 obj.axes_plot{ii}
                util.plot.inner_title(['clip' num2str(ii)], 'position', 'bottom', 'ax', obj.axes_plot{ii});
                
            end
            
        end
        
        function update(obj,~,~)
            
            if ~obj.checkGUI
                return;
            end
            
            %%%%%%%%%%%%%%%%%% panel controls %%%%%%%%%%%%%%%%%%%%%%
            
            obj.input_cut_size.String = ['cut size= ' num2str(obj.clip.cut_size)];
            obj.input_cut_size.FontSize = obj.edit_font_size;
            
            if obj.clip.use_mex
                obj.button_mex_cut.String = 'mex cut';
            else
                obj.button_mex_cut.String = 'no mex cut';
            end
            
            obj.button_mex_cut.FontSize = obj.font_size;
            
            %%%%%%%%%%%%%%%%%% panel adjust %%%%%%%%%%%%%%%%%%%%%%%%
            
            if obj.clip.use_lock_adjust
                obj.button_use_lock_adjust.String = 'adjust lock';
            else                
                obj.button_use_lock_adjust.String = 'no lock';
            end
                        
            obj.button_use_lock_adjust.FontSize = obj.font_size;
            
            if obj.clip.use_moments
                obj.button_use_moments.String = 'use moments';
            else
                obj.button_use_moments.String = 'brightest pix';
            end
            
            obj.button_use_moments.FontSize = obj.font_size;
                        
            %%%%%%%%%%%%%%%%%% panel find stars %%%%%%%%%%%%%%%%%%%%
            
            obj.input_num_stars.String = ['num stars= ' num2str(obj.clip.num_stars)];
            obj.input_num_stars.FontSize = obj.edit_font_size;
            
            obj.input_filter_size.String = ['filter size= ' num2str(obj.clip.filter_size)];
            obj.input_filter_size.FontSize = obj.edit_font_size;
            
            obj.input_filter_sigma.String = ['filter sigma= ' num2str(obj.clip.filter_sigma)];
            obj.input_filter_sigma.FontSize = obj.edit_font_size;
            
            obj.input_avoid_edges.String = ['avoid edges= ' num2str(obj.clip.avoid_edges)];
            obj.input_avoid_edges.FontSize = obj.edit_font_size;
            
            %%%%%%%%%%%%%%%%%% panel plot %%%%%%%%%%%%%%%%%%%%%%%%%%
            
%             plot(obj.axes_plot, (1:size(obj.clip.mean_drift,1))', obj.clip.mean_drift);
                        
            %%%%%%%%%%%%%%%%%% panel positions %%%%%%%%%%%%%%%%%%%%%
            
            str = '';
            
            for ii = 1:size(obj.clip.positions,1)
                
                str = [str 'x= ' num2str(obj.clip.positions(ii,1))];
                str = [str ' | y= ' num2str(obj.clip.positions(ii,2)) char(10)];
                
            end
            
            obj.button_positions.String = str;
            obj.button_positions.FontSize = obj.edit_font_size;
            
            %%%%%%%%%%%%%%%%%% panel close %%%%%%%%%%%%%%%%%%%%%%%%%
            
            obj.button_close.FontSize = obj.font_size;
            
        end
        
        function updateGUI(obj)
            obj.update;
        end
        
        function c = check(obj)
            
            c = obj.checkGUI;
            
        end
        
        function c = checkGUI(obj)
           
            c = ~isempty(obj.panel_plot) && isvalid(obj.panel_plot);
            
        end
        
    end
    
    methods % callbacks
        
        %%%%%%%%%%%%%%%%%% panel controls %%%%%%%%%%%%%%%%%%%%%%

        function callback_cut_size(obj, hndl, ~)
                      
            numbers = util.text.extract_numbers(hndl.String);
            
            numbers = cell2mat(numbers);
            
            if isempty(numbers)
                numbers = obj.clip.default_cut_size;
            end
            
            if obj.debug_bit>1, disp(['callback: cut_size= ' num2str(numbers)]); end
            
            obj.clip.cut_size = numbers;
            
            obj.update;
            
        end
        
        function callback_mex_cut(obj, ~, ~)
           
            if obj.debug_bit>1, disp('callback: mex_cut'); end
            
            obj.clip.use_mex = ~obj.clip.use_mex;
            
            obj.update;
            
        end
        
        %%%%%%%%%%%%%%%%%% panel adjust %%%%%%%%%%%%%%%%%%%%%%%%

        function callback_use_lock_adjust(obj, ~, ~)
           
            if obj.debug_bit>1, disp('callback: use_lock_adjust'); end
            
            obj.clip.use_lock_adjust = ~obj.clip.use_lock_adjust;
            
            obj.update;
            
        end
        
        function callback_use_moments(obj, ~, ~)
           
            if obj.debug_bit>1, disp('callback: use_moments'); end
            
            obj.clip.use_moments = ~obj.clip.use_moments;
            
            obj.update;
            
        end
        
        %%%%%%%%%%%%%%%%%% panel find stars %%%%%%%%%%%%%%%%%%%%

        function callback_num_stars(obj, hndl, ~)
                      
            numbers = util.text.extract_numbers(hndl.String);
            
            numbers = cell2mat(numbers);
            
            if isempty(numbers)
                numbers = 1;
            end
            
            if obj.debug_bit>1, disp(['callback: num_stars= ' num2str(numbers)]); end
            
            obj.clip.num_stars = numbers;
            
            obj.update;
            
        end
        
        function callback_filter_size(obj, hndl, ~)
                      
            numbers = util.text.extract_numbers(hndl.String);
            
            numbers = cell2mat(numbers);
            
            if isempty(numbers)
                numbers = obj.clip.default_filter_size;
            end
            
            if obj.debug_bit>1, disp(['callback: filter_size= ' num2str(numbers)]); end
            
            obj.clip.filter_size = numbers;
            
            obj.update;
            
        end
        
        function callback_filter_sigma(obj, hndl, ~)
                      
            numbers = util.text.extract_numbers(hndl.String);
            
            numbers = cell2mat(numbers);
            
            if isempty(numbers)
                numbers = obj.clip.default_filter_sigma;
            end
            
            if obj.debug_bit>1, disp(['callback: filter_sigma= ' num2str(numbers)]); end
            
            obj.clip.filter_sigma = numbers;
            
            obj.update;
            
        end
        
        function callback_avoid_edges(obj, hndl, ~)
                      
            numbers = util.text.extract_numbers(hndl.String);
            
            numbers = cell2mat(numbers);
            
            if isempty(numbers)
                numbers = ceil(obj.clip.cut_size/2);
            end
            
            if obj.debug_bit>1, disp(['callback: avoid_edges= ' num2str(numbers)]); end
            
            obj.clip.avoid_edges = numbers;
            
            obj.update;
            
        end
        
        %%%%%%%%%%%%%%%%%% panel close %%%%%%%%%%%%%%%%%%%%%%%%%

        function callback_close(obj, ~, ~)
           
            if obj.debug_bit>1, disp('callback: close'); end
            
            delete(obj.fig.fig);
            
        end
        
    end
    
    
    
end