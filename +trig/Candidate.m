classdef Candidate < handle

    properties(Transient=true)
        
    end
    
    properties % objects
        
        head; % header information
        kernel_props; % struct with the data about the best kernel
        star_props; % table with one row from astrometry/GAIA
        
    end
    
    properties % inputs/outputs
        
        identification_time = ''; % timestring when this event was identified by analysis (UTC)
        
        notes = {}; % anything added to the event, like failed cuts
        
        cutouts; % cutouts of the selected star, for the duration of the extended region
        stack; % stack of the images in which this event is detected (optional)
        
        time_index; % out of the extended region
        kern_index; % from the full filter bank
        star_index; % from the stars that passed the initial burn-in (not from the subset that survived the pre-filter
        star_index_global; % from all the stars in this run (compare this to the catalog, cutouts, etc)
        frame_index; % the frame of the peak inside the original batch/file
                
        time_range; % time indices around time_index that are considered "part of the event"
        
        filenames; % cell array with the filenames for each frame in the extended region
        frame_numbers; % frame numbers for each frame inside its respective file
        
        % all of these are for the duration of the extended region
        timestamps; 
        juldates;
        flux_raw; % all fluxes for all stars, without any processing
        flux_corrected; % fluxes for all stars, corrected by PSD or by removing linear fit
        flux_filtered; % flux after matched-filtering, for the peak star and kernel
        
        auxiliary; % values of auxiliary measurements (e.g., background, offsets) for all stars and all frames in the extended region
        aux_names; % cell array with the names of each auxiliary
        aux_indices; % struct with a field for each auxiliary name, containing the index in the big matrix
        
        used_psd_corr; % did we apply a PSD correction?
        used_filt_std; % did we correct for each filtered flux's noise after filtering?
        
        snr;
        threshold;
        
    end
    
    properties % switches/controls
        
        kept = 1; % by default any new event is considered real (kept)
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        
        
    end
    
    properties(Hidden=true)
        
        is_simulated = 0; % by default events are not simulated
        sim_pars; % struct with the simulation parameters (for simulated events only)

        % lower thresholds used to find time range, kernels and stars around
        % the peak that are also included. 
        thresh_time;
        thresh_kern;
        thresh_star;
        
        kern_extra; % any other kernels that passed the lower threshold for kernels
        star_extra; % any other stars that passed the lower threshold for stars
        
        version = 1.00;
        
    end
    
    methods % constructor
        
        function obj = Candidate(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'trig.Candidate')
                if obj.debug_bit>1, fprintf('Candidate copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit>1, fprintf('Candidate constructor v%4.2f\n', obj.version); end
            
            end
            
        end
        
    end
    
    methods % reset/clear
        
    end
    
    methods % getters
        
    end
    
    methods % setters
        
        function addNote(obj, str)
            
            if isempty(obj.notes)
                obj.notes{1} = str;
            else
                obj.notes{end+1} = str;
            end
            
        end
        
    end
    
    methods % calculations
        
    end
    
    methods % plotting tools / GUI
        
    end    
    
end

