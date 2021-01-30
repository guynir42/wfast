classdef StarMonitor < handle
% Collect some additional info about specific stars during an analysis run.
% This object is part of trig.EventFinder, and will always get the photometric
% data + processed cutouts for each batch. If the star_indices prop is
% empty then it does nothing. 
% If there are some star indices then it will save the flux and auxiliary
% data for that star for all batches in the run. 
% If any events trigger at this star it will also save those Candidate
% objects but with added "cutouts_all" that includes cutouts for all stars
% in the time period of the event. 
%
% NOTE: the data is only collected after the burn-in period. The star
%       indices are therefore in the reduced list of stars. Make sure the
%       index you give is from that list, not from the list of all stars!

    properties(Transient=true)
        
    end
    
    properties % objects
        
        cand@trig.Candidate;
        
    end
    
    properties % inputs/outputs
        
        timestamps; % timestamp of each measurement
        juldates; % julian date of each measurement
        flux; % flux for each star, for the entire run
        aux; % additional photometric properties like offsets and width and background
        
    end
    
    properties % switches/controls
        
        star_indices = []; % list the indices of stars you want to monitor, out of the stars passing the threshold after the burn-in
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        
        
    end
    
    properties(Hidden=true)
       
        version = 1.00;
        
    end
    
    methods % constructor
        
        function obj = StarMonitor(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'trig.StarMonitor')
                if obj.debug_bit>1, fprintf('StarMonitor copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit>1, fprintf('StarMonitor constructor v%4.2f\n', obj.version); end
            
            end
            
        end
        
    end
    
    methods % reset/clear
        
        function reset(obj)
            
            obj.cand = trig.Candidate.empty;
            
            obj.timestamps = [];
            obj.juldates = []; 
            obj.flux = []; 
            obj.aux = []; 
            
        end
        
    end
    
    methods % getters
        
    end
    
    methods % setters
        
    end
    
    methods % calculations
        
        function input(obj, timestamps, juldates, flux, aux, cutouts, candidates)
            
            if isempty(obj.star_indices)
                return;
            end
            
            
            obj.timestamps = vertcat(obj.timestamps, timestamps); 
            obj.juldates = vertcat(obj.juldates, juldates); 
            obj.flux = vertcat(obj.flux, flux(:,obj.star_indices)); 
            obj.aux = vertcat(obj.aux, aux(:,obj.star_indices,:)); 
            
            for ii = 1:length(obj.star_indices)
                for jj = 1:length(candidates)
                    if candidates(jj).star_index==obj.star_indices(ii) % note: this index is on the stars that passed the threshold after the burn in
                        c = util.oop.full_copy(candidates(jj));
                        c.cutouts_all = cutouts;
                        obj.cand = vertcat(obj.cand, c);
                    end
                end
            end
            
        end
        
    end
    
    methods % plotting tools / GUI
        
    end    
    
end

