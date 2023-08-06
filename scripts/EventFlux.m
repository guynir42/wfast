classdef EventFlux < handle

    properties
        
        filename; 
        base_dir = fullfile(getenv('DATA'), 'WFAST\two_year_results\occultation_flux');
        raw_data;
        latest_plot_type = '';
        latest_figure = [];
        
    end
    
    properties(Dependent=true)
        cand;
        mcmc;
        flux;
        time;
        frame_index;
        star_index;
        aperture_index;
        event_date;
        
    end
    
    methods % constructor
        
        function obj = EventFlux(event_date, base_dir)
            
            if nargin < 1 || isempty(event_date)
                event_date = [];
            end
            
            if nargin < 2 || isempty(base_dir)
                folder = obj.base_dir;
            else
                folder = base_dir;
                obj.base_dir = base_dir;
            end
            
            if isempty(event_date)
                obj.filename = uigetfile({'*.mat'}, 'Select a file', folder);
            else
                obj.filename = sprintf('occult_%s_fluxes.mat', event_date); 
            end
            
        end
        
    end
    
    methods % getters
        
        function val = get.cand(obj)
            if isempty(obj.raw_data)
                val = [];
            else
                val = obj.raw_data.occultation;
            end
        end
        
        function val = get.mcmc(obj)
            if isempty(obj.raw_data) || isempty(obj.raw_data.occultation)
                val = [];
            else
                val = obj.raw_data.occultation.mcmc;
            end
        end
        
        function val = get.flux(obj)
            if isempty(obj.raw_data)
                val = [];
            else
                val = obj.raw_data.flux;
            end
        end
        
        function val = get.time(obj)
            if isempty(obj.raw_data)
                val = [];
            else
                val = obj.raw_data.time;
            end
        end
        
        function val = get.frame_index(obj)
            if isempty(obj.raw_data)
                val = [];
            else
                val = obj.raw_data.frame_index;
            end
        end
        
        function val = get.star_index(obj)
            if isempty(obj.raw_data)
                val = [];
            else
                val = obj.raw_data.star_index;
            end
        end
        
        function val = get.aperture_index(obj)
            if isempty(obj.raw_data)
                val = [];
            else
                val = obj.raw_data.aperture_index;
            end
        end
        
        function val = get.event_date(obj)
            if isempty(obj.filename)
                val = [];
            else
                idx = regexp(obj.filename, '\d{4}-\d{2}-\d{2}');
                if isempty(idx)
                    val = [];
                else
                    val = obj.filename(idx:idx+9); 
                end
            end
        end
        
        function val = getVelocity(obj)
            
            if isempty(obj.cand) || isempty(obj.head)
                val = [];
            else
                val = obj.cand.head.ephem.getShadowVelocity(1); 
            end
            
        end
        
    end
    
    methods % loading/saving/calculations
        
        function load(obj)
            
            if isempty(obj.filename)
                error('Cannot load data without a filename!')
            end
            
            obj.raw_data = load(fullfile(obj.base_dir, obj.filename)); 
            
        end
        
        function save(obj)
            
            occultation = obj.cand;
            flux = obj.flux;
            time = obj.time;
            frame_index = obj.frame_index;
            star_index = obj.star_index;
            aperture_index = obj.aperture_index;
            save(fullfile(obj.base_dir, obj.filename), 'occultation', 'flux', 'time', 'frame_index', 'star_index', 'aperture_index', '-v7.3')
            
        end
        
        function print(obj)
            
            if isempty(obj.latest_plot_type) || isempty(obj.latest_figure) || ~isvalid(obj.latest_figure)
                error('Cannot print figure, no plot type or figure handle given'); 
            end
            
            name = sprintf('%s_%s', obj.latest_plot_type, obj.event_date);
            fig = gcf;
            if ~strcmp(fig.Name, name)
                error('Current figure name ("%s") does not match filename "%s"', fig.Name, name);
            end
            
            util.sys.print(fullfile(fileparts(obj.base_dir), 'plots', name));
            
        end
        
    end
    
    methods % plotting tools
        
        function fig = showFluxCutouts(obj)
            
            obj.latest_plot_type = 'flux_cutouts'; 
            
            fh = util.plot.FigHandler(sprintf('%s_%s', obj.latest_plot_type, obj.event_date));
            fh.width = 30;
            fh.height = 15;
            fh.clear;

            obj.cand.showFluxCutouts('parent', fh.fig);
            
            if nargout > 0
                fig = fh.fig;
            end
            
            obj.latest_figure = fh.fig;
            
        end
        
        function fig = showNeighbors(obj, number, add_offsets)
            
            if nargin<2 || isempty(number)
                number = [];
            end
            
            if nargin<3 || isempty(add_offsets)
                add_offsets = false; 
            end
            
            obj.latest_plot_type = 'neighbors'; 
            
            fh = util.plot.FigHandler(sprintf('%s_%s', obj.latest_plot_type, obj.event_date));
            fh.width = 30;
            fh.height = 15;
            fh.clear;

            obj.cand.showNearestStars('parent', fh.fig, 'number', number, 'offsets', add_offsets);
            
            if nargout > 0
                fig = fh.fig;
            end
            
            obj.latest_figure = fh.fig;
            
        end
        
        function fig = showMCMC(obj)
            
            obj.latest_plot_type = 'mcmc'; 
            
            fh = util.plot.FigHandler(sprintf('%s_%s', obj.latest_plot_type, obj.event_date));
            fh.width = 30;
            fh.height = 15;
            fh.clear;

            obj.mcmc.showResults('parent', fh.fig); 
            
            if nargout > 0
                fig = fh.fig;
            end
            
            obj.latest_figure = fh.fig;
            
        end
        
        function fig = showOutliers(obj)
            
            obj.latest_plot_type = 'outliers'; 
            
            fh = util.plot.FigHandler(sprintf('%s_%s', obj.latest_plot_type, obj.event_date));
            fh.width = 30;
            fh.height = 15;
            fh.clear;
            
            idx1 = obj.frame_index - 2000;
            if idx1 < 1, edx1 = 1; end
            
            idx2 = obj.frame_index + 2000;
            if idx2 > size(obj.flux, 1), idx2 = size(obj.flux, 1); end
            
            f = obj.flux(idx1:idx2, :, obj.aperture_index);
            obj.cand.showOutlierAnalysis(f, 'parent', fh.fig, 'recalc', 1); 
            
            if nargout > 0
                fig = fh.fig;
            end
            
            obj.latest_figure = fh.fig;
            
        end
        
    end
    
end



