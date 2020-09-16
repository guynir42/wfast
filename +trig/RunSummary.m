classdef RunSummary < handle
% Lightweight object that summarizes the parameters, the data quality and 
% the key results from an occultation scanning run. 
% It is produced by the EventFinder using the produceSummary() method. 
% 
% The parameters for each relevant object are saved as structs, "finder_pars", 
% "store_pars", and "checker_pars". The header info is save in "head" as usual. 
% 
% An important result is the histogram of star-time per S/N bin. 
% This is stored without separating into bins for each star, but only for 
% each S/N value. This is done so we can stack the total hours or bin them
% by e.g., ecliptic latitude later on. 
% 
% Other interesting results are the "snr_values" saved for each batch in the 
% run, the cut values (again, histogrammed by values only, not by star), and 
% the parameters of failed/passed simulated events. 
% That last one can be used to find what occultation parameters manage to 
% pass our threshold (and it can be compared to the theoretical values). 
% accumulating these events can help define the real thresholds for different
% types of occultations. 
% 
    
    properties % objects
        
        head;
        
        finder_pars; 
        store_pars;
        checker_pars;
        
    end
    
    properties % inputs/outputs
        
        snr_values; 
        runtime; 
        total_batches; 
        
        snr_bin_edges; % the edges used in the star_seconds and losses histograms
        star_seconds; % the number of useful seconds accumulated each S/N bin (this is saved after subtracting losses)
        star_seconds_with_losses; % the number of seconds accumulated without excluding anything
        losses_exclusive; % same as histogram, only counting the losses due to each cut (exclusive means times where ONLY this cut was responsible)
        losses_inclusive; % same as histogram, only counting the losses due to each cut (inclusive means times where this cut was ALSO responsible)
        losses_bad_stars; % number of seconds lost when disqualifying stars from the black list 
        
        cut_names = {}; % get these from the QualityChecker
        cut_indices = {}; % get these from the QualityChecker
        cut_thresholds = []; % get these from the QualityChecker
        cut_two_sided = []; % get these from the QualityChecker
        
        cut_histograms; % accumualted values of different cuts for each value, summed over all stars (dim1 is values, dim2 is cut type)
        cut_bin_edges; % a vector with the bin edges for each cut (same dimensions as the histograms)
        
        black_list_stars; % a vector of star numbers (from the full list of stars?)
        good_stars; % a list of the stars that passed the threshold after the burn-in period
        
        % structs describing simulated events that passed/failed the threshold
        sim_events;
        
    end
    
    properties % switches/controls
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        
        
    end
    
    properties(Hidden=true)
        
        version = 1.00;
        
    end
    
    properties(Hidden=true, Transient=true)
        
        results; 
        axes_names;
        axes_indices; 
        type_result=''; 
        R_hist;
        R_edges;
        r_hist;
        r_edges;
        b_hist;
        b_edges;
        v_hist;
        v_edges;
        
        
    end
    
    methods % constructor
        
        function obj = RunSummary(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'trig.RunSummary')
                if obj.debug_bit>1, fprintf('RunSummary copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit>1, fprintf('RunSummary constructor v%4.2f\n', obj.version); end
            
            end
            
        end
        
    end
    
    methods % reset/clear
        
    end
    
    methods % getters
        
    end
    
    methods % setters
        
    end
    
    methods % calculations
       
        function val = calcStats(obj, varargin)
            
            import util.text.cs; 
            
            input = util.text.InputVars;
            input.input_var('result', 'fraction'); % can choose "fraction" or "snr"
            input.input_var('star_size_bin', 0.1); % size of bin for R (star radius)
            input.input_var('radius_bin', 0.1); % size of bin for r (occulter radius)
            input.input_var('impact_bin', 0.2); % size of bin for b (impact parameter)
            input.input_var('velocity_bin', 5); % size of bin for v (velocity)
            input.input_var('x_axis', 'r'); % which parameter to put in the x-axis of each subframe (dim 2)
            input.input_var('y_axis', 'b'); % which parameter to put in the y-axis of each subframe (dim 1)
            input.input_var('horizontal', 'v'); % which parameter to put in the horizontal spread of frames (dim 3)
            input.input_var('vertical', 'R'); % which parameter to put in the vertical spread of frames (dim 4)
            input.scan_vars(varargin{:}); 
            
            val = ''; 
            
            %%%%%%%%%%% find the bin edges and distributions of each variable %%%%%%%%%%%%%%%%%
            if cs(input.result, 'fraction')
                events = obj.sim_events; % mix the two together (the failed events have detect_snr = []); 
                mode = 1; 
            elseif cs(input.result, 'snr')
                events = obj.sim_events(logical([obj.sim_events.passed])); % take only the passing events for S/N calculation 
                mode = 2; 
            else
                error('Unkown option "%s" to input "result". Try using "fraction" or "snr". ', input.result); 
            end
            
            [obj.R_hist, obj.R_edges] = histcounts([events.R]+0.001, 'BinWidth', input.star_size_bin);
            [obj.r_hist, obj.r_edges] = histcounts([events.r]+0.001, 'BinWidth', input.radius_bin);
            [obj.b_hist, obj.b_edges] = histcounts([events.b]+0.001, 'BinWidth', input.impact_bin);
            [obj.v_hist, obj.v_edges] = histcounts([events.v]+0.001, 'BinWidth', input.velocity_bin);

            %%%%%%%%%% setup the size of the matrix for the result %%%%%%%%%%%%%%%%%
            
            % name of the dimension (single letter)
            dim1 = input.y_axis;
            dim2 = input.x_axis;
            dim3 = input.horizontal;
            dim4 = input.vertical;
            
            obj.axes_names = {dim1,dim2,dim3,dim4}; 
            obj.axes_indices.dim1 = dim1;
            obj.axes_indices.dim2 = dim2;
            obj.axes_indices.dim3 = dim3;
            obj.axes_indices.dim4 = dim4;
            
            % edge vectors
            e1 = obj.([dim1 '_edges']); 
            e2 = obj.([dim2 '_edges']); 
            e3 = obj.([dim3 '_edges']); 
            e4 = obj.([dim4 '_edges']); 
            
            % add one more bin 
%             e1 = [e1 2*e1(end)-e1(end-1)];
%             e2 = [e2 2*e2(end)-e2(end-1)];
%             e3 = [e3 2*e3(end)-e3(end-1)];
%             e4 = [e4 2*e4(end)-e4(end-1)];
            
            % adjust the bins just so they don't fall on the uniformly sampled grid:
            e1 = e1 - 0.001;
            e2 = e2 - 0.001;
            e3 = e3 - 0.001;
            e4 = e4 - 0.001;
            
            % number of bins (equal to number of edges minus one)
            Ndim1 = length(e1)-1; 
            Ndim2 = length(e2)-1; 
            Ndim3 = length(e3)-1; 
            Ndim4 = length(e4)-1; 
            
            val = zeros(Ndim1, Ndim2, Ndim3, Ndim4, 'single'); 
            
            %%%%%%%%%% start binning the data %%%%%%%%%%%%%%%%%%

            for ii = 1:Ndim1
                
                subset1 = events([events.(dim1)]>=e1(ii) & [events.(dim1)]<e1(ii+1)); % choose only events inside this range
                
                for jj = 1:Ndim2
                    
                    subset2 = subset1([subset1.(dim2)]>=e2(jj) & [subset1.(dim2)]<e2(jj+1)); % choose only events inside this range
                    
                    for kk = 1:Ndim3
                        
                        subset3 = subset2([subset2.(dim3)]>=e3(kk) & [subset2.(dim3)]<e3(kk+1)); % choose only events inside this range
                    
                        for mm = 1:Ndim4
                            
                            subset4 = subset3([subset3.(dim4)]>=e4(mm) & [subset3.(dim4)]<e4(mm+1)); % choose only events inside this range

                            if mode==1 % fraction mode
                                
                                if isempty(subset4)
                                    val(ii,jj,kk,mm) = NaN; % don't have a better idea what to put here
                                else
                                    passed = nnz([subset4.passed]); 
                                    total = numel(subset4); 

                                    val(ii,jj,kk,mm) = passed./total; 
                                end
                                
                            elseif mode==2 % snr mode
                                
                                if isempty(subset4)
                                    val(ii,jj,kk,mm) = 0; % don't have a better idea what to put here
                                else
                                    val(ii,jj,kk,mm) = mean([subset4.detect_snr]./[subset4.star_snr].*10); % adjust all S/N values to a standard star with S/N=10
                                end
                                
                            end % mode

                        end % for mm
                    
                    end % for kk
                    
                end % for jj
                
            end % for ii
            
            obj.results = val; 
            obj.type_result = input.result; % keep track of what we calculated
            
            if nargin==0
                clear val;
            end
                                    
        end
        
    end
    
    methods % plotting tools / GUI
        
        function showSimulations(obj, varargin)
           
            input = util.text.InputVars;
%             input.input_var('result', 'fraction'); % can choose "fraction" or "snr"
%             input.input_var('star_size_bin', 0.1); % size of bin for R (star radius)
%             input.input_var('radius_bin', 0.1); % size of bin for r (occulter radius)
%             input.input_var('impact_bin', 0.2); % size of bin for b (impact parameter)
%             input.input_var('velocity_bin', 5); % size of bin for v (velocity)
%             input.input_var('x_axis', 'r'); % which parameter to put in the x-axis of each subframe (dim 2)
%             input.input_var('y_axis', 'b'); % which parameter to put in the y-axis of each subframe (dim 1)
%             input.input_var('horizonatal', 'v'); % which parameter to put in the horizontal spread of frames (dim 3)
%             input.input_var('vertical', 'R'); % which parameter to put in the vertical spread of frames (dim 4)
            input.input_var('parent', []); % which axes to plot to? default is gca()
            input.input_var('font_size', 18); % fonts on the axes
            input.scan_vars(varargin{:}); 
            
            if isempty(input.parent)
                input.parent = gcf;
            end
            
            delete(input.parent.Children);
            
            obj.calcStats(varargin{:}); % pass through all the relevant inputs
            
            % get shorthands for the edges and names of the variables
            name1 = obj.axes_names{1}; % b
            edges1 = obj.([name1 '_edges']); % b edges
            
            name2 = obj.axes_names{2}; % r
            edges2 = obj.([name2 '_edges']); % r edges
            
            name3 = obj.axes_names{3}; % v
            edges3 = obj.([name3 '_edges']); % v edges
            
            name4 = obj.axes_names{4}; % R
            edges4 = obj.([name4 '_edges']); % R edges
            
            % parameters for the axes grid
            gap = 0.01; 
            margin_x = 0.08;
            margin_y = 0.20;
            Nx = length(edges3)-1;
            Ny = length(edges4)-1;
            width =  (1 - margin_x)/Nx - gap;
            height = (1 - margin_y)/Ny - gap;
            
            ax = {}; 
            
            for ii = 1:Ny % run over R values
                
                for jj = 1:Nx % run over v values
                    
                    ax{ii,jj} = axes('Parent', input.parent, ...
                        'Position', [margin_x + (jj-1)*(width+gap) margin_y + (ii-1)*(height+gap) width height]); 
                    
                    if strcmp(obj.type_result, 'fraction')
                        util.plot.show(obj.results(:,:,jj,ii), 'monochrome', 1,'fancy', 'on', ...
                            'x_values', edges2(1:end-1), 'y_values', edges1(1:end-1)); 
                        
                        colorbar(ax{ii,jj}, 'off'); 
                        axis(ax{ii,jj}, 'normal'); 
                        ax{ii,jj}.YDir = 'normal';
                    elseif strcmp(obj.type_result, 'snr')
                        
                        contour(ax{ii,jj}, edges2(1:end-1),edges1(1:end-1), obj.results(:,:,jj,ii), [7.5 10 12.5 15 17.5 20], 'k', 'ShowText', 'on'); 
                        
                    end
                    
                    if ii==1
                        
                    else
                        
                    end
                    
                    if ii==Ny
                        [n,u] = obj.getFullName(name2);
                        xlabel(ax{ii,jj}, sprintf('%s \n[%s]', n,u));
                        ax{ii,jj}.XTick = [edges2(2) edges2(floor(length(edges2)/2)) edges2(end-3)]; 
                    else
                        ax{ii,jj}.XTick = []; 
                        
                        if ~isempty(ax{ii,jj}.YTick)
                            ax{ii,jj}.YTick(end) = [];
                        end
                        
                    end
                    
                    if jj==1
                        [n,u] = obj.getFullName(name1);
                        ylabel(ax{ii,jj}, sprintf('%s [%s]', n,u));
                    else
                        ax{ii,jj}.YTick = [];
                    end
                    
                    if jj==Nx
                        
                    else
                        if ~isempty(ax{ii,jj}.XTick)
%                             ax{ii,jj}.XTick(end) = [];
                        end
                    end
                    
                    
                    ax{ii,jj}.FontSize = input.font_size; 
                    
                    util.plot.inner_title(sprintf('%s= %4.2f \n %s= %4.2f', name4, edges4(ii), name3, edges3(jj)), ...
                        'Position', 'NorthWest', 'FontSize', input.font_size, 'margin', 0.06); 
                    
                end % for jj (v)
                
            end % for ii (R)
                        
        end
        
        function [name, units] = getFullName(obj, letter)
            
            if letter=='R'
                name = 'stellar size'; 
            elseif letter=='r'
                name = 'radius';
            elseif letter=='b'
                name = 'impact parameter'; 
            elseif letter=='v'
                name = 'velocity';
            end
            
            if nargout>1
                if letter=='R'
                    units = 'FSU'; 
                elseif letter=='r'
                    units = 'FSU';
                elseif letter=='b'
                    units = 'FSU'; 
                elseif letter=='v'
                    units = 'FSU/s';
                end 
            end
            
        end
        
    end    
    
end

