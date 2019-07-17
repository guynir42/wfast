classdef FilterBank < handle

    properties(Transient=true)
        
    end
    
    properties % objects
        
        gen@occult.CurveGenerator;
        prog@util.sys.ProgressBar;
        
    end
    
    properties % inputs/outputs
        
        % inputs: 
        R_list = 0:0.25:0.5; % star radius, FSU
        r_list = 0.1:0.1:2; % occulter radius, FSU
        b_list = 0:0.2:2; % impact parameter, FSU
        v_list = 3:1:30; % crossing velocity (projected), FSU/second
        % lets assume t=0 for all!
        
        W = 4; % time window, seconds
        T = 30; % integration time, ms
        f = 25; % frame rate, Hz
        
        % outputs:
        timestamps; % a 1D time axis for all LCs
        bank; % a 5D map, with lightcurves in 1st dimension, and parameters in the rest... 
        snr_map; % a 4D map, for each parameter combination, what is the minimal S/N with all its neighbors. 
        
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
        
        function obj = FilterBank(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'occult.FilterBank')
                if obj.debug_bit, fprintf('FilterBank copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            
            elseif ~isempty(varargin) && isa(varargin{1}, 'occult.CurveGenerator')
                if obj.debug_bit, fprintf('FilterBank generator-based constructor v%4.2f\n', obj.version); end
                
                obj.gen = varargin{1};
            
                obj.prog = util.sys.ProgressBar;
                
            else
                if obj.debug_bit, fprintf('FilterBank constructor v%4.2f\n', obj.version); end
                
                obj.gen = occult.CurveGenerator;
                
                obj.prog = util.sys.ProgressBar;
                
            end
            
        end
        
    end
    
    methods % reset/clear
        
        function reset(obj)
            
            obj.gen.reset;
            
            obj.bank = [];
            obj.timestamps = [];
            obj.snr_map = [];
            
        end
        
    end
    
    methods % getters
        
        function val = num_pars(obj)
            
            val = length(obj.R_list).*length(obj.r_list).*length(obj.b_list).*length(obj.v_list);
            
        end
        
        function val = mem_size_gb(obj)
            
            val = obj.num_pars.*length(obj.timestamps)*8/1024^3; 
            
        end
        
        function val = getNeighborsList(obj)
            
            val = zeros(0,4);
            
            for ii = 1:4
                
                for jj = 1:2
                    
                    val(end+1,ii) = (-1).^jj;
                    
                end
                
            end
            
        end
        
    end
    
    methods % setters
        
    end
    
    methods % calculations
        
        function makeBank(obj)
            
            obj.bank = [];
            
            obj.prog.start(length(obj.R_list));
            
            % because loops in matlab are fun! 
            for ii = 1:length(obj.R_list)
                
                for jj = 1:length(obj.r_list)
                    
                    for kk = 1:length(obj.b_list)
                        
                        for mm = 1:length(obj.v_list)
                            
                            obj.gen.R = obj.R_list(ii);
                            obj.gen.r = obj.r_list(jj);
                            obj.gen.b = obj.b_list(kk);
                            obj.gen.v = obj.v_list(mm);
                            
                            obj.gen.getLightCurves;
                            obj.bank(:,ii,jj,kk,mm) = obj.gen.lc.flux;
                            
                        end % for mm (v_list)
                        
                    end % for kk (b_list)
                    
                end % for jj (r_list)
                
                obj.prog.show(ii);
                
            end % for ii (R_list)
            
            obj.timestamps = obj.gen.lc.time;
            
        end
        
        function checkBank(obj)
            
            obj.prog.start(length(obj.R_list));
            
            obj.snr_map = [];
            
            NL = obj.getNeighborsList;
            
            for ii = 1:length(obj.R_list)
                
                for jj = 1:length(obj.r_list)
                    
                    for kk = 1:length(obj.b_list)
                        
                        for mm = 1:length(obj.v_list)
                            
                            values = [];
                            
                            for nn = 1:size(NL,1) % 4 dimensional neighborhood
                                
%                                 idx = [ii+(-1).^(nn) jj+(-1).^floor(nn/2) kk+(-1).^floor(nn/4) mm+(-1).^floor(nn/8)];

                                idx = NL(nn,:) + [ii,jj,kk,mm];
                                
                                if any(idx<1)
                                    continue;
                                end
                                
                                if idx(1)>length(obj.R_list) || idx(2)>length(obj.r_list) || idx(3)>length(obj.b_list) || idx(4)>length(obj.v_list)
                                    continue;
                                end
                                
                                values = [values obj.compareLightcurves(obj.bank(:,ii,jj,kk,mm), obj.bank(:,idx(1), idx(2), idx(3), idx(4)))];
                                
                            end
                            
                            if isempty(values)
                                obj.snr_map(ii,jj,kk,mm) = NaN;
                            else
                                obj.snr_map(ii,jj,kk,mm) = max(values);
                            end
                            
                        end % for mm (v_list)
                        
                    end % for kk (b_list)
                    
                end % for jj (r_list)
                
                obj.prog.show(ii);
                
            end % for ii (R_list) 
            
        end
        
        function snr = compareLightcurves(obj, this_flux, that_flux)
            
            f1 = this_flux - 1;
            f2 = that_flux - 1;
            
            % now the normalized kernels
            k1 = f1./sqrt(sum(f1.^2));
            k2 = f2./sqrt(sum(f2.^2));
            
            self_signal = sum(k1.*f1);
            cross_signal = sum(k2.*f1);
            
            snr = cross_signal./self_signal;
            
        end
        
    end
    
    methods % plotting tools / GUI
        
        function plotLC(obj, idx)
            
            ax = gca;
            
            cla(ax);
            
            core_lc = obj.bank(:, idx(1), idx(2), idx(3), idx(4));
            
            h = plot(ax, obj.timestamps, core_lc);
            h.DisplayName = sprintf('R(%d) | r(%d) | b(%d) | v(%d)', idx(1), idx(2), idx(3), idx(4));
            
            hold(ax, 'on');
            
            NL = obj.getNeighborsList;
            
            for nn = 1:size(NL,1) % 4 dimensional neighborhood
                                
%                 idx = [R_idx+(-1).^(nn) r_idx+(-1).^floor(nn/2) b_idx+(-1).^floor(nn/4) v_idx+(-1).^floor(nn/8)];
                neigh_idx = idx + NL(nn,:);

                if any(neigh_idx<1)
                    continue;
                end

                if neigh_idx(1)>length(obj.R_list) || neigh_idx(2)>length(obj.r_list) || neigh_idx(3)>length(obj.b_list) || neigh_idx(4)>length(obj.v_list)
                    continue;
                end

                new_lc = obj.bank(:, neigh_idx(1), neigh_idx(2), neigh_idx(3), neigh_idx(4));
                loss = obj.compareLightcurves(core_lc, new_lc);
                
                h = plot(ax, obj.timestamps, new_lc, ':');
                h.DisplayName = sprintf('R(%d) | r(%d) | b(%d) | v(%d) | loss= %f', neigh_idx(1), neigh_idx(2), neigh_idx(3), neigh_idx(4), loss);
                
            end
            
            hold(ax, 'off');
            
            hl = legend(ax, 'Location', 'best');
            
            hl.FontSize = 10;
            
        end
        
    end    
    
end

