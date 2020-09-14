classdef EvFinderGUI < handle
    
    properties 
        
        owner@trig.EventFinder; % link back to containg object

        fig@util.plot.FigHandler;
             
        buttons = {};
        menus = {};
        panels = {};
        
        latest_error = '';
        latest_warning = '';
        
        font_size = 11;
        big_font_size = 14;
        edit_font_size = 10;
        small_font_size = 9;
        
        color_on = [0 0.3 1];
        
        debug_bit = 1;
        
    end
    
    properties % gui stuff
        
        menu_pars;
        
        panel_controls;
        panel_summary;
        panel_plots;
        panel_sim; 
        
        
        panel_close;
        button_close;
        
        panel_image;
        button_reset_axes;
        axes_image;
    
    end
    
    properties (Hidden=true)
              
        version = 1.00;
        
    end
            
    methods % constructor
       
        function obj = EvFinderGUI(owner)
            
            % later add other options like copy constructor
            if isa(owner, 'trig.EventFinder')
                
                if obj.debug_bit>1, fprintf('EvFinderGUI constructor v%4.2f\n', obj.version); end
                
                obj.owner = owner;
                
            else
                error('Input an EventFinder to constructor of EvFinderGUI!');
            end
            
        end
        
    end
    
    methods % make/update/check GUI
    
        function make(obj)
            
            import util.plot.GraphicButton;
            import util.plot.GraphicPanel;
            import util.plot.ContrastLimits;
            import util.plot.MenuItem;
            
            obj.buttons = {};
            obj.menus = {};
            obj.panels = {}; 
            
            obj.fig = util.plot.FigHandler('EventFinder');
            obj.fig.clear;
            obj.fig.bottom = 5;
            obj.fig.height = 20;
            obj.fig.width = 40;
            obj.fig.center;
%             obj.fig.maximize;
            
            
            %%%%%%%%%%%%%%%%%%%%%%% MENUS %%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % MenuItem(parent, text, type, variable, tooltip, separator)
            % obj.addButton(name, text, type, variable, tooltip, separator)
            % menu types: menu, toggle, push, input, input_text, info, custom
            
            
            %%%%%%%%%%% pars menu %%%%%%%%%%%%%%%
            
            obj.menu_pars = MenuItem(obj, '&Parameters', 'menu'); 
            
            obj.menu_pars.addButton('menu_thresh', '&Thresholds', 'menu'); 
            obj.menu_pars.menu_thresh.addButton('input_thresh', '&Main threshold', 'input', 'pars.threshold', 'main event detection threshold');
            obj.menu_pars.menu_thresh.addButton('input_time', '&Time threshold', 'input', 'pars.time_range_thresh', 'threshold for continuous time range around peak (absolute or relative value)', 1);
            obj.menu_pars.menu_thresh.addButton('input_kern', '&Kern threshold', 'input', 'pars.kern_range_thresh', 'threshold for other kernels around peak (absolute or relative value)');
            obj.menu_pars.menu_thresh.addButton('input_star', '&Star threshold', 'input', 'pars.star_range_thresh', 'threshold for other stars around peak (absolute or relative value)');
            obj.menu_pars.menu_thresh.addButton('input_spread', '&Minimal spread', 'input', 'pars.min_time_spread', 'minimal number of frames on either side of peak that count as part of the event', 1);
            
            obj.menu_pars.addButton('menu_num_events', '&Number events', 'menu'); 
            obj.menu_pars.menu_num_events.addButton('input_max_events', '&Max events', 'input', 'pars.max_events', 'maximum number of event finding iterations on each batch'); 
            obj.menu_pars.menu_num_events.addButton('input_limit_batch', '&Batch limit', 'input', 'pars.limit_events_per_batch', 'number of events for a single batch that would put it on the black list'); 
            obj.menu_pars.menu_num_events.addButton('input_limit_star', '&Star limit', 'input', 'pars.limit_events_per_star', 'number of events for a single star that would put it on the black list'); 
            
            obj.menu_pars.addButton('menu_corrections', '&Corrections', 'menu'); 
            obj.menu_pars.menu_corrections.addButton('button_cr', 'remove &CRs', 'toggle', 'store.pars.use_remove_cosmic_rays', 'remove cosmic ray spikes from the flux measurements'); 
            obj.menu_pars.menu_corrections.addButton('input_cr', 'CR &Threshold', 'input', 'store.pars.cosmic_ray_threshold', 'threshold for removing cosmic rays'); 
            obj.menu_pars.menu_corrections.addButton('button_psd', 'use &PSD', 'toggle', 'pars.use_psd_correction', 'apply Power Spectral Density correction to the fluxes', 1); 
            obj.menu_pars.menu_corrections.addButton('button_std', 'use filtered &STD', 'toggle', 'pars.use_std_filtered', 'correct the filtered flux by the STD calculated on the background region'); 
            obj.menu_pars.menu_corrections.addButton('button_keep_var', 'keep &Variances', 'toggle', 'pars.use_keep_variances', 'keep the values of the filtered fluxes variance for all batches', 1); 
            
            obj.menu_pars.addButton('menu_filters', '&Filters', 'menu'); 
            obj.menu_pars.menu_filters.addButton('button_prefilter', '&Prefilter', 'toggle', 'pars.use_prefilter', 'filter on a small template-bank with lower threshold and only use the full bank on passing stars');
            obj.menu_pars.menu_filters.addButton('input_prefilter', '&Threshold', 'input', 'pars.pre_threshold', 'threshold for the smaller template bank');
            obj.menu_pars.menu_filters.addButton('input_bank', '&Filter bank', 'input_text', 'pars.filter_bank_full_filename', 'path and name of the full filter bank MAT-file, relative to DATA', 1);
            obj.menu_pars.menu_filters.addButton('input_small_bank', '&Small bank', 'input_text', 'pars.filter_bank_small_filename', 'path and name of the small filter bank MAT-file, relative to DATA');
            
            obj.menu_pars.addButton('menu_simulations', '&Simulations', 'menu');
            obj.menu_pars.menu_simulations.addButton('button_sim', 'use &Simulations', 'toggle', 'pars.use_sim', 'turn on/off the injection simulations'); 
            obj.menu_pars.menu_simulations.addButton('input_num_sim', '&Num events', 'input', 'pars.num_sim_events_per_batch', 'how many simulated events (on average) for each batch (can be fractional)');
            
            %%%%%%%%%%% store menu %%%%%%%%%%%%%%%
            
            obj.menu_pars.addButton('menu_store', '&DataStore', 'menu', '', '', 1); 
            obj.menu_pars.menu_store.addButton('input_length_psd', 'length &PSD', 'input', 'store.pars.length_psd', 'maximum buffer size (frames) used for calculating the PSD'); 
            obj.menu_pars.menu_store.addButton('input_length_bg', 'length B/&G', 'input', 'store.pars.length_background', 'number of frames of the background buffer used for calculating the rms'); 
            obj.menu_pars.menu_store.addButton('input_length_ext', 'length &Extended', 'input', 'store.pars.length_extended', 'number of frames of the extended region used for filtering overlap'); 
            obj.menu_pars.menu_store.addButton('input_length_search', 'length &Search', 'input', 'store.pars.length_search', 'number of frames of the search region used for finding events');
            obj.menu_pars.menu_store.addButton('input_length_burn_in', 'length &Burn in', 'input', 'store.pars.length_burn_in', 'number of frames used in the burn-in period before starting the event search');
            
            obj.menu_pars.menu_store.addButton('button_threshold', '&Use star threshold', 'toggle', 'store.pars.use_threshold', 'remove events with low S/N after finishing the burn-in', 1);
            obj.menu_pars.menu_store.addButton('input_threshold', 'star &Threshold', 'input', 'store.pars.threshold', 'threshold for removing events');
            
            obj.menu_pars.menu_store.addButton('input_aperture', 'Aperture radius', 'info', 'store.aperture_radius', '', 1);
            
            %%%%%%%%%%% checker menu %%%%%%%%%%%%%%%
            
            obj.menu_pars.addButton('menu_checker', '&QualityChecker', 'menu'); 
            
            obj.menu_pars.menu_checker.addButton('menu_apply', '&Apply cuts', 'menu'); 
            obj.menu_pars.menu_checker.menu_apply.addButton('button_delta_t', 'delta &T', 'toggle', 'store.checker.pars.use_delta_t', 'disqualify based on irregular jumps in the timestamps'); 
            obj.menu_pars.menu_checker.menu_apply.addButton('button_shakes', '&Shakes', 'toggle', 'store.checker.pars.use_shakes', 'disqualify based on flux-weighted mean offset size'); 
            obj.menu_pars.menu_checker.menu_apply.addButton('button_defocus', '&Defocus', 'toggle', 'store.checker.pars.use_defocus', 'disqualify based on the flux-weighted mean PSF width'); 
            obj.menu_pars.menu_checker.menu_apply.addButton('button_slope', '&Flux slope', 'toggle', 'store.checker.pars.use_slope', 'disqualify based on global linear changes in the average flux'); 
            obj.menu_pars.menu_checker.menu_apply.addButton('button_bad_rows', '&Bad rows/cols', 'toggle', 'store.checker.pars.use_near_bad_rows_cols', 'disqualify based on proximity to bad rows/columns'); 
            obj.menu_pars.menu_checker.menu_apply.addButton('button_offset_size', '&Offset size', 'toggle', 'store.checker.pars.use_offset_size', 'disqualify based on large offset size'); 
            obj.menu_pars.menu_checker.menu_apply.addButton('button_linear_motion', '&Linear motion', 'toggle', 'store.checker.pars.use_linear_motion', 'disqualify based on linear changes in the x/y offsets'); 
            obj.menu_pars.menu_checker.menu_apply.addButton('button_bg', '&Background', 'toggle', 'store.checker.pars.use_background_intensity', 'disqualify based on background count'); 
            obj.menu_pars.menu_checker.menu_apply.addButton('button_nan_flux', '&NaN flux', 'toggle', 'store.checker.pars.use_nan_flux', 'disqualify based on NaN values in the flux', 1); 
            obj.menu_pars.menu_checker.menu_apply.addButton('button_nan_offsets', 'NaN &XY', 'toggle', 'store.checker.pars.use_nan_offsets', 'disqualify based on NaN values in the x/y offsets'); 
            obj.menu_pars.menu_checker.menu_apply.addButton('button_photo_flag', '&Photo-flag', 'toggle', 'store.checker.pars.use_photo_flag', 'disqualify based on any flag raised by the photometry pipeline');            
            obj.menu_pars.menu_checker.menu_apply.addButton('button_corr', '&Correlations', 'toggle', 'store.checker.pars.use_correlations', 'disqualify based on correlations of the flux to some auxiliary data', 1); 
            
            obj.menu_pars.menu_checker.addButton('menu_thresh', '&Thresholds', 'menu'); 
            obj.menu_pars.menu_checker.menu_thresh.addButton('input_delta_t', 'delta &T', 'input', 'store.checker.pars.thresh_delta_t', 'threshold for irregular jumps in the timestamps'); 
            obj.menu_pars.menu_checker.menu_thresh.addButton('input_shakes', '&Shakes', 'input', 'store.checker.pars.thresh_shakes', 'threshold for flux-weighted mean offset size'); 
            obj.menu_pars.menu_checker.menu_thresh.addButton('input_defocus', '&Defocus', 'input', 'store.checker.pars.thresh_defocus', 'threshold for the flux-weighted mean PSF width'); 
            obj.menu_pars.menu_checker.menu_thresh.addButton('input_slope', '&Flux slope', 'input', 'store.checker.pars.thresh_slope', 'threshold for global linear changes in the average flux'); 
            obj.menu_pars.menu_checker.menu_thresh.addButton('input_offset_size', '&Offset size', 'input', 'store.checker.pars.thresh_offset_size', 'threshold for large offset size'); 
            obj.menu_pars.menu_checker.menu_thresh.addButton('input_linear_motion', '&Linear motion', 'input', 'store.checker.pars.thresh_linear_motion', 'threshold for linear changes in the x/y offsets'); 
            obj.menu_pars.menu_checker.menu_thresh.addButton('input_bg', '&Background', 'input', 'store.checker.pars.thresh_background_intensity', 'threshold for background count');          
            obj.menu_pars.menu_checker.menu_thresh.addButton('input_corr', '&Correlations', 'input', 'store.checker.pars.thresh_correlation', 'threshold for correlations of the flux to some auxiliary data', 1); 
            
            obj.menu_pars.menu_checker.addButton('button_dilate', '&Dilate', 'toggle', 'store.checker.pars.use_dilate', 'expand the time around each bad frame', 1); 
            obj.menu_pars.menu_checker.addButton('input_dilate', 'dilate &Region', 'input', 'store.checker.pars.dilate_region', 'expand by this number of frames in either direction'); 
            
            obj.menu_pars.menu_checker.addButton('button_subtract_offsets', 'subtract &Offsets', 'toggle', 'store.checker.pars.use_subtract_mean_offsets', 'subtract the star-average offset X/Y from each individual star offset', 1); 
            obj.menu_pars.menu_checker.addButton('input_smoothing_slope', '&Smoothing slope', 'input', 'store.checker.pars.smoothing_slope', 'time-scale for slope cut (frames)'); 
            obj.menu_pars.menu_checker.addButton('input_dist_bad', '&Distance bad rows/cols', 'input', 'store.checker.pars.distance_bad_rows_cols', 'events closer than this many pixels to a bad row/column are disqualified');
            obj.menu_pars.menu_checker.addButton('input_linear_timescale', '&Linear timescale', 'input', 'store.checker.pars.linear_timescale', 'time-scale for linear motion cut (frames)'); 
            
            obj.menu_pars.menu_checker.addButton('input_num_edges', '&Num hist edges', 'input', 'store.checker.pars.num_hist_edges', 'how many edges/bins to span the cut values histograms', 1); 
            
            %%%%%%%%%%% hours menu %%%%%%%%%%%%%%%
            
            obj.menu_pars.addButton('menu_hours', 'Star &Hours', 'menu'); 
            obj.menu_pars.menu_hours.addButton('input_bin_min', '&Low S/N edge', 'input', 'store.checker.hours.snr_bin_min', 'low edge of S/N histogram');
            obj.menu_pars.menu_hours.addButton('input_bin_max', '&High S/N edge', 'input', 'store.checker.hours.snr_bin_max', 'high edge of S/N histogram');
            obj.menu_pars.menu_hours.addButton('input_bin_width', 'S/N &Width', 'input', 'store.checker.hours.snr_bin_width', 'edge width of S/N histogram');
            
            obj.menu_pars.assignJavaObjectsTopLevel; 
            
            %%%%%%%%%%%%%%%%%%% LEFT SIDE %%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            N_left = 13; % number of buttons on left side
            
            pos = N_left;
            
            %%%%%%%%%%% panel controls %%%%%%%%%%%%%%%
            
            % Add buttons using obj.addButton(button_name, var_name='', type='', str1='', str2='', font_size='', split=1, color_on=[], color_off=[], tooltip)
            
            num_buttons = 3;
            pos = pos-num_buttons;
            obj.panel_controls = GraphicPanel(obj.owner, [0 pos/N_left 0.2 num_buttons/N_left], 'controls', 1); % last input is for vertical (default)
            obj.panel_controls.number = num_buttons;
            obj.panel_controls.addButton('input_threshold', 'pars.threshold', 'input', 'thresh= ', '', '', 1, [], [], 'detection threshold for the main filter bank'); 
            obj.panel_controls.addButton('button_psd', 'pars.use_psd_correction', 'toggle', 'PSD off', 'PSD on', '', 0.5, obj.color_on, 'red', 'use Power Spectra Density correction of flux'); 
            obj.panel_controls.addButton('button_std', 'pars.use_std_filtered', 'toggle', 'filt. std off', 'filt std on', '', 0.5, obj.color_on, 'red', 'use filtered flux STD correction'); 
            obj.panel_controls.addButton('button_cr_remove', 'store.pars.use_remove_cosmic_rays', 'toggle', 'CR remove off', 'CR remove on', '', 0.5, obj.color_on, 'red', 'remove cosmic rays spikes'); 
            obj.panel_controls.addButton('input_cr_thresh', 'store.pars.cosmic_ray_threshold', 'input', 'CR thresh= ', '', '', 0.5, [], [], 'threshold for removal of cosmic rays'); 
            obj.panel_controls.margin = [0.02 0.02]; 
            obj.panel_controls.make;
            
            
            %%%%%%%%%%% panel summary %%%%%%%%%%%%%%%
            
            num_buttons = 4;
            pos = pos-num_buttons;
            obj.panel_summary = GraphicPanel(obj.owner, [0 pos/N_left 0.2 num_buttons/N_left], 'summary', 1); % last input is for vertical (default)
            obj.panel_summary.number = num_buttons;
            obj.panel_summary.addButton('button_runtime', 'runtime_batches_str', 'info', ' ', '', '', 1, [], [], 'total search time and number of batches processed (after burn-in)'); 
            obj.panel_summary.addButton('button_star_hours', 'star_hours_str', 'info', ' ', '', '', 1, [], [], 'useful/total star hours'); 
            obj.panel_summary.addButton('button_num_stars', 'num_stars_str', 'info', ' ', '', '', 1, [], [], 'number of stars passing the threshold out of total number of cutouts'); 
            obj.panel_summary.addButton('button_cand', 'num_candidates', 'info', 'cand= ', '', '', 0.5, [], [], 'number of triggered event candidates'); 
            obj.panel_summary.addButton('button_kept', 'num_kept', 'info', 'kept= ', '', '', 0.5, [], [], 'number of candidates that passed all cuts');             
            obj.panel_summary.margin = [0.02 0.02]; 
            obj.panel_summary.make;
            
            %%%%%%%%%%% panel plots %%%%%%%%%%%%%%%
            
            num_buttons = 4;
            pos = pos-num_buttons;
            obj.panel_plots = GraphicPanel(obj.owner, [0 pos/N_left 0.2 num_buttons/N_left], 'plots', 1); % last input is for vertical (default)
            obj.panel_plots.number = num_buttons;
            obj.panel_plots.addButton('button_cand', 'showCandidates', 'push', 'cand', '', '', 0.5, [], [], 're-render the candidate viewer'); 
            obj.panel_plots.addButton('button_psd', 'popupPSD', 'push', 'PSD', '', '', 0.5, [], [], 'pop up a window with the PSD plot'); 
            obj.panel_plots.addButton('button_stars', 'popupStars', 'push', 'stars', '', '', 0.5, [], [], 'pop up a window with the star S/N plot'); 
            obj.panel_plots.addButton('button_snr', 'popupSNR', 'push', 'S/N', '', '', 0.5, [], [], 'pop up a window with a histogram of the best S/N measured in eacj batch'); 
            obj.panel_plots.addButton('button_quality', 'popupRunQuality', 'push', 'quality', '', '', 0.5, [], [], 'pop up a window with a plot of the airmass, background and PSF width over the entire run');
            obj.panel_plots.addButton('button_cuts', 'popupCuts', 'push', 'cuts', '', '', 0.5, [], [], 'pop up a window with the quality cuts viewer'); 
            obj.panel_plots.addButton('button_sim', 'popupSimulated', 'push', 'sim', '', '', 0.5, [], [], 'pop up a window with the results of injected events processing'); 
            obj.panel_plots.addButton('button_star_hours', 'popupStarHours', 'push', 'star hours', '', '', 0.5, [], [], 'pop up a window with the star hour viewer'); 
            
            obj.panel_plots.margin = [0.02 0.02]; 
            obj.panel_plots.make;
            
            
            %%%%%%%%%%% panel sim %%%%%%%%%%%%%%%
            
            num_buttons = 1;
            pos = pos-num_buttons;
            obj.panel_sim = GraphicPanel(obj.owner, [0 pos/N_left 0.2 num_buttons/N_left], 'simulations', 1); % last input is for vertical (default)
            obj.panel_sim.number = num_buttons;
            obj.panel_sim.addButton('button_use_sim', 'pars.use_sim', 'toggle', 'sim off', 'sim on', '', 0.5, obj.color_on, 'red', 'turn on/off the injection simulations'); 
            obj.panel_sim.addButton('input_num_sim', 'pars.num_sim_events_per_batch', 'input', 'N= ', '', '', 0.5, [], [], 'number of events to be injected to each batch, on average (can be fractional!)'); 
            
            obj.panel_sim.margin = [0.02 0.02]; 
            obj.panel_sim.make;
            
            
            
            %%%%%%%%%%% panel image %%%%%%%%%%%%%%%%%%
            
            obj.panel_image = uipanel('Title', '', 'Position', [0.2 0 0.8 1]);
            
            %%%%%%%%%%% panel close %%%%%%%%%%%%%%%%%%
            
            obj.panel_close = uipanel('Position', [0 0 0.2 1/N_left]);
            obj.button_close = GraphicButton(obj.panel_close, [0 0 1 1], obj.owner, '', 'custom', 'CLOSE GUI');
            obj.button_close.Callback = @obj.callback_close;
            
            obj.owner.showCandidates;
            
            obj.update;
            
        end
        
        function update(obj,~,~)
                        
            if ~obj.check
                return;
            end
           
            for ii = 1:length(obj.buttons)
                obj.buttons{ii}.update;
            end
            
            for ii = 1:length(obj.menus)
                obj.menus{ii}.update;
            end
            
        end
                        
        function c = check(obj)
           
            c = ~isempty(obj) && ~isempty(obj.panel_close) && isvalid(obj.panel_close);
            
        end
        
    end
                
    methods % callbacks
        
        function callback_close(obj, ~, ~)
           
            if obj.debug_bit>1, disp('callback: close'); end
            
            delete(obj.fig.fig);
            
        end
        
    end
    
end