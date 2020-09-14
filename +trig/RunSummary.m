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
    
    properties(Transient=true)
        
    end
    
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
        sim_events_passed;
        sim_events_failed;
        
    end
    
    properties % switches/controls
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        
        
    end
    
    properties(Hidden=true)
       
        version = 1.00;
        
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
        
    end
    
    methods % plotting tools / GUI
        
    end    
    
end

