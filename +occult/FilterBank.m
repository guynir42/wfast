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
        r_list = 0.3:0.1:2; % occulter radius, FSU
        b_list = 0:0.2:2; % impact parameter, FSU
        v_list = 5:5:30; % crossing velocity (projected), FSU/second
        % lets assume t=0 for all!
        
        W = 4; % time window, seconds
        T = 30; % integration time, ms
        f = 25; % frame rate, Hz
        
        % outputs:
        timestamps; % a 1D time axis for all LCs
        bank; % a 5D map, with lightcurves in 1st dimension, and parameters in the rest... 
        snr_map; % a 4D map, for each parameter combination, what is the minimal S/N with all its neighbors. 
        
        sim_snr; % simulation result S/N after dividing by star S/N (5D matrix with zero where there was no detection)
        filtered_bank; % all templates, after filtering them with all kernels for event finding
        filtered_index; % auxiliary index to tell where we are in the last 4 dimensions of "filtered_bank"
        
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
            
            obj.snr_map = [];
            obj.sim_snr = [];
            
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
        
        function val = sim_pars(obj)
            
            if isempty(obj.filtered_index) || obj.filtered_index==0
                val = [];
            else
                
                S = size(obj.bank);
                S = S(2:end);
                [R,r,b,v] = ind2sub(S,obj.filtered_index);

                val = struct('R', obj.R_list(R), 'r', obj.r_list(r), 'b', obj.b_list(b), 'v', obj.v_list(v), 'W', obj.W, 'T', obj.T, 'f', obj.f); 

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
                            obj.bank(:,ii,jj,kk,mm) = single(obj.gen.lc.flux);
                            
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
                                
                                values = [values obj.compareKernels(obj.bank(:,ii,jj,kk,mm), obj.bank(:,idx(1), idx(2), idx(3), idx(4)))];
                                
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
        
        function snr = compareKernels(obj, this_flux, that_flux, full_xcorr)
            
            if nargin<4 || isempty(full_xcorr)
                full_xcorr = 0;
            end
            
            f1 = this_flux - 1;
            f2 = that_flux - 1;
            
            % now the normalized kernels
            k1 = f1./sqrt(sum(f1.^2));
            k2 = f2./sqrt(sum(f2.^2));
            
            if full_xcorr==0
                self_signal = sum(k1.*f1);
                cross_signal = sum(k2.*f1);
            else
                self_signal = max(filter2(f1,k1));
                cross_signal = max(filter2(f1,k2));
            end
            
            snr = cross_signal./self_signal;
            
        end
        
        function values = monteCarloCheck(obj, num_trials, full_xcorr)
           
            if nargin<3 || isempty(full_xcorr)
                full_xcorr = 1;
            end
            
            obj.prog.start(num_trials);
            
            values = NaN(num_trials,1);
            
            flat_bank = reshape(obj.bank, [size(obj.bank,1), obj.num_pars]);
            
            for ii = 1:num_trials
                
                obj.gen.lc.pars = obj.randomPars;
                obj.gen.getLightCurves;
                
                snr = util.stat.maxnd(occult.compareKernels(obj.gen.lc.flux, flat_bank, full_xcorr));
                
                if snr<0.9
                    fprintf('S/N = %f | R= %f | r= %f | b= %f | v= %f\n', snr, obj.gen.R, obj.gen.r, obj.gen.b, obj.gen.v);
                end
                
                values(ii) = snr;
                
                obj.prog.showif(ii);
                
            end
            
        end
        
        function pars = randomPars(obj)
            
            pars = occult.Parameters;
            
            vec = obj.R_list;
            pars.R = rand*(max(vec)-min(vec))+min(vec);
            
            vec = obj.r_list;
            pars.r = rand*(max(vec)-min(vec))+min(vec);
            
            vec = obj.b_list;
            pars.b = rand*(max(vec)-min(vec))+min(vec);
            
            vec = obj.v_list;
            pars.v = rand*(max(vec)-min(vec))+min(vec);
            
        end
        
        function score = calcMinStarSNR(obj, varargin)
            
            if isempty(obj.bank)
                error('Must fill the filter bank first!');
            end
            
            obj.prog.start(length(obj.R_list));
            
            % because loops in matlab are fun! 
            for ii = 1:length(obj.R_list)
                
                for jj = 1:length(obj.r_list)
                    
                    for kk = 1:length(obj.b_list)
                        
                        for mm = 1:length(obj.v_list)
                            
                        end % for mm (v_list)
                        
                    end % for kk (b_list)
                    
                end % for jj (r_list)
                
                obj.prog.show(ii);
                
            end % for ii (R_list)
            
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
                loss = occult.compareKernels(core_lc, new_lc);
                
                h = plot(ax, obj.timestamps, new_lc, ':');
                h.DisplayName = sprintf('R(%d) | r(%d) | b(%d) | v(%d) | loss= %f', neigh_idx(1), neigh_idx(2), neigh_idx(3), neigh_idx(4), loss);
                
            end
            
            hold(ax, 'off');
            
            hl = legend(ax, 'Location', 'best');
            
            hl.FontSize = 10;
            
        end
        
    end    
    
end

