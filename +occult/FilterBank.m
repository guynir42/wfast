classdef FilterBank < handle
% Contains a set of filter kernels, each corresponding to one set of 
% occultation parameters, that can be used for matched-filtering or for 
% parameter estimation. 
% The parameters are sorted in 4 dimensions: R, r, b and v (set t=0). 
% The range and step for these parameters are defined in "R_list", "r_list", 
% "b_list" and "v_list". 
% Additional parameters are t=0, W, T and f that are shared by all kernels. 
%
% The "bank" property actually contains all the lightcurves. 
% This is a 5D matrix where the 1st dimension is the time axis, and the 
% rest are for each of the occultation parameters. 
%
% This object can calculate the analytical S/N for each kernel using the 
% findAnalyticalSNR(), checkNeighbors() and checkMonteCarlo(). 
% 
% This class is also used to calculate S/N from the trig.Finder object's 
% simulations, injecting each kernel into real data and storing the trigger
% result in "snr_sim_full". 
% Plotting these results is done using showSNR(). 


    properties(Transient=true)
        
    end
    
    properties % objects
        
        gen@occult.CurveGenerator; % this produces the kernel bank
        prog@util.sys.ProgressBar; % track the time it takes to calculate things
        
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
        
        filtered_bank; % all templates, after filtering them with all kernels for event finding
        filtered_index; % auxiliary index to tell where we are in the last 4 dimensions of "filtered_bank"
        
        snr_neighbors; % a 4D map, for each parameter combination, what is the minimal S/N with all its neighbors. 
        snr_self_test; % a 1D list of S/N results for self test
        snr_analytical; % a 4D map for each kernel what is the detection S/N for a signal with the same parameters 
        snr_sim_full; % simulation result S/N after dividing by star S/N (2D matrix with zero where there was no detection)
        
    end
    
    properties(Dependent=true)
        
        snr_simulated; % a 4D map, for each parameter combination, the S/N for detection on real data (the median of "snr_sim_full" reshaped to 4D). 

    end
    
    properties % switches/controls
        
        threshold = 0.95; % minimal overlap to be considered close enough to another kernel
        
        debug_bit = 1;
        
    end
    
    properties(Hidden=true)
       
        version = 1.01;
        
    end
    
    methods % constructor
        
        function obj = FilterBank(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'occult.FilterBank')
                if obj.debug_bit>1, fprintf('FilterBank copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            
            elseif ~isempty(varargin) && isa(varargin{1}, 'occult.CurveGenerator')
                if obj.debug_bit>1, fprintf('FilterBank generator-based constructor v%4.2f\n', obj.version); end
                
                obj.gen = varargin{1};
            
                obj.prog = util.sys.ProgressBar;
                
            else
                if obj.debug_bit>1, fprintf('FilterBank constructor v%4.2f\n', obj.version); end
                
                obj.gen = occult.CurveGenerator;
                
                obj.prog = util.sys.ProgressBar;
                
            end
            
        end
        
    end
    
    methods % reset/clear
        
        function reset(obj)
            
            obj.snr_analytical = [];
            obj.snr_sim_full = [];
            
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
        
        function val = get.snr_simulated(obj)
           
            if isempty(obj.snr_sim_full)
                val = [];
            else
                
                val = median(obj.snr_sim_full,1,'omitnan');

                val = reshape(val, [length(obj.R_list), length(obj.r_list), length(obj.b_list), length(obj.v_list)]);
                
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
        
        function checkNeighbors(obj)
            
            obj.prog.start(length(obj.R_list));
            
            obj.snr_analytical = [];
            
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
                                
                                values = [values obj.compareKernels(single(obj.bank(:,ii,jj,kk,mm)-1), single(obj.bank(:,idx(1), idx(2), idx(3), idx(4))-1))];
                                
                            end
                            
                            if isempty(values)
                                obj.snr_neighbors(ii,jj,kk,mm) = NaN;
                            else
                                obj.snr_neighbors(ii,jj,kk,mm) = max(values);
                            end
                            
                        end % for mm (v_list)
                        
                    end % for kk (b_list)
                    
                end % for jj (r_list)
                
                obj.prog.show(ii);
                
            end % for ii (R_list) 
            
        end
        
        function checkMonteCarlo(obj, num_trials, full_xcorr)
           
            if nargin<3 || isempty(full_xcorr)
                full_xcorr = 1;
            end
            
            obj.prog.start(num_trials);
            
            obj.snr_self_test = NaN(num_trials,1);
            
            flat_bank = reshape(obj.bank, [size(obj.bank,1), obj.num_pars]);
            
            for ii = 1:num_trials
                
                obj.gen.lc.pars = obj.randomPars;
                obj.gen.getLightCurves;
                
                snr = max(occult.compareKernels(single(flat_bank-1), single(obj.gen.lc.flux-1), full_xcorr));
                
                if snr<obj.threshold
                    fprintf('S/N = %f | R= %f | r= %f | b= %f | v= %f\n', snr, obj.gen.R, obj.gen.r, obj.gen.b, obj.gen.v);
                end
                
                obj.snr_self_test(ii) = snr;
                
                obj.prog.showif(ii);
                
            end
            
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
        
        function [lc, pars] = randomLC(obj, fix_R)
            
            if nargin<2 || isempty(fix_R)
                fix_R = [];
            end
            
            S = size(obj.bank); 
            
            if isempty(fix_R)
                
                idx = randperm(prod(S(2:end)),1); 
                
                [sub_R, sub_r, sub_b, sub_v] = ind2sub(S(2:end), idx);
                
                lc = obj.bank(:,idx); 
                    
            else
                
                idx = randperm(prod(S(3:end)),1);
                
                sub_R = find(fix_R==obj.R_list);
                
                if isempty(sub_R)
                    error('Wrong value of fix_R. Choose from the R_list= %s', util.text.print_vec(obj.R_list)); 
                end
                
                [sub_r, sub_b, sub_v] = ind2sub(S(3:end), idx);
                
                lc = obj.bank(:, sub_R, sub_r, sub_b, sub_v); 
                
                idx = idx + prod(S(3:end))*sub_R; 
                
            end
                
            pars.R = obj.R_list(sub_R); 
            pars.r = obj.r_list(sub_r); 
            pars.b = obj.b_list(sub_b); 
            pars.v = obj.v_list(sub_v); 
            pars.W = obj.W;
            pars.T = obj.T;
            pars.f = obj.f; 
            
            pars.bank_index = idx; 
            
        end
        
        function score = calcMinStarSNR(obj, varargin) % need to finish this! 
            
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
        
        function findAnalyticalSNR(obj, kernels)
            
            obj.snr_analytical = [];
            
            f = single(obj.bank-1); % "raw" lightcurves (not including noise)
                
            if nargin<2 || isempty(kernels)
                
                k = f./sqrt(sum(f.^2)); % kernels as they would be used by the filtering code

                obj.snr_analytical = max(util.vec.convolution(k,f), 'conj', 1); 
            
            else
                
                k = kernels; 
                k = k./sqrt(sum(kernels.^2));
                
                % because loops in matlab are fun! 
                for ii = 1:length(obj.R_list)

                    for jj = 1:length(obj.r_list)

                        for kk = 1:length(obj.b_list)

                            for mm = 1:length(obj.v_list)

                                obj.snr_analytical(1,ii,jj,kk,mm) = util.stat.maxnd(util.vec.convolution(k,f(:,ii,jj,kk,mm), 'conj', 1)); 
                                
                            end % for mm (v_list)

                        end % for kk (b_list)

                    end % for jj (r_list)

                end % for ii (R_list)

            end
                
            obj.snr_analytical = permute(obj.snr_analytical, [2,3,4,5,1]);
            
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
        
        function showSNR(obj, varargin)
            
            import util.text.cs;
            
            input = util.text.InputVars;
            input.input_var('parent', []);
            input.input_var('data', 'analytical'); % can also choose "simulated"
            input.input_var('mode', 'contour'); % can also choose "heat"
            input.input_var('star', 10, 'star_snr');
            input.input_var('intervals', 10:5:40);
            input.input_var('threshold', 7.5);
            input.input_var('font_size', 20);
            input.scan_vars(varargin{:});
            
            if isnumeric(input.data)
                snr = input.data;
            elseif cs(input.data, 'analytical')
                snr = obj.snr_analytical;
            elseif cs(input.data, 'simulated')
                snr = obj.snr_simulated;
            else
                error('Input "data" as numeric S/N matrix or choose "analytical" or "simulated" options'); 
            end
                
            if isempty(snr)
                error('S/N data is empty!');
            end
            
            if isempty(input.parent)
                input.parent = gcf;
            end
            
            delete(input.parent.Children);
            
            snr = permute(snr, [3,2,1,4]); % now snr has order: b,r,R,v (for better plotting)

            NR = length(obj.R_list);
            Nv = length(obj.v_list);
            
            margin = 0.02;
            
            x_start = 0.07;
            x_step = (1-x_start)/Nv;
            
            y_start = 0.15;
            y_step = (1-y_start)/NR;
            
            ax = {};
            counter = 1;

            ypos = y_start;

            for ii = 1:NR

                xpos = x_start;

                for jj = 1:Nv

                    ax{counter} = axes('Parent', input.parent, 'Position', [xpos ypos, x_step-margin, y_step-margin]);

                    if cs(input.mode, 'contour')
                        
                        contour(obj.r_list, obj.b_list, snr(:,:,ii,jj)*input.star, input.intervals, 'k', 'ShowText', 'on');
                        
                        ax{counter}.NextPlot = 'add';
                        
                        if input.threshold>0
                            contour(obj.r_list, obj.b_list, snr(:,:,ii,jj)*input.star, input.threshold.*[1 1], 'k--', 'ShowText', 'on');
                        end
                        
                        ax{counter}.NextPlot = 'replace';
                        
                    elseif cs(input.mode, 'heatmap')
                        util.plot.show(snr(:,:,ii,jj)*input.star, 'xvalues', obj.r_list, 'yvalues', obj.b_list, 'fancy', 'off');
                    else
                        error('Unknown mode option "%s". Choose "contour" or "heatmap"', input.mode);
                    end

                    ax{counter}.XTick = [0.5 1 1.5 2];
                    
                    vec = ax{counter}.XTick;
                    
                    for k = 1:length(vec)
                        if ii==1 && vec(k)<2
                            ax{counter}.XTickLabels{k} = num2str(vec(k));
                        else
                            ax{counter}.XTickLabels{k} = '';
                        end

                        ax{counter}.XAxis.MinorTick = 'on';

                        grid(ax{counter}, 'on');

                    end
                    
                    if ii==1, xlabel(ax{counter}, 'radius'); end

                    ax{counter}.YTick = [0 0.5 1 1.5 2];
                    
                    vec = ax{counter}.YTick;
                    
                    for k = 1:length(vec)
                        if jj==1 && (ii==1 || vec(k)>0) % && ( ii==NR || vec(k)<2 ) && ( ii==1 || vec(k)>0 )
                            ax{counter}.YTickLabels{k} = num2str(vec(k));
                        else
                            ax{counter}.YTickLabels{k} = '';
                        end

                        ax{counter}.YAxis.MinorTick = 'on';

                        grid(ax{counter}, 'on');

                    end
                    
                    if jj==1
                        ylabel(ax{counter}, 'impact par.');
                    end

                    if jj==Nv && cs(input.mode, 'heatmap') 
                        
                        colorbar(ax{counter}, 'on');
                        drawnow;
                        ax{counter}.Position(3) = ax{counter-1}.Position(3);
                        ax{counter}.Position(4) = ax{counter-1}.Position(4);
                        
                    end
                    
                    ax{counter}.FontSize = input.font_size;
                    h = util.plot.inner_title(ax{counter}, sprintf('R= %4.2f | v= %d', obj.R_list(ii), obj.v_list(jj)), ...
                        'Position', 'NorthWest', 'Color', 'black', 'FontSize', input.font_size-4); 

                    counter = counter + 1;

                    xpos = xpos + x_step;

                end
                
                ypos = ypos + y_step;

            end

            
        end
        
    end    
    
end

