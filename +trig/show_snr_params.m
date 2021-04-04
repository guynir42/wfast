function h = show_snr_params(ev,varargin)
% Usage: h = show_snr_params(ev,varargin)

    if nargin==0, help('trig.show_snr_params'); return; end
    
    if ~isa(ev, 'struct') || ~isfield(ev, 'D') || ~isfield(ev, 'R') || ...
            ~isfield(ev,'r') || ~isfield(ev, 'b') || ~isfield(ev, 'v') || ...
            ~isfield(ev, 'detect_snr')
        error('Must input a struct array of simulated events with fields: "D", "R", "r", "b", "v", "passed", and "detect_snr"'); 
    end
    
    input = util.text.InputVars;
    input.input_var('pars', 'rb', 'par_order', 'parameter_order'); 
    input.input_var('distance', 40); 
    input.input_var('axes', [], 'axis'); 
    input.input_var('marker_size', 20); 
    input.input_var('font_size', 18); 
    input.scan_vars(varargin{:}); 
    
    if isempty(input.axes)
        input.axes = gca;
    end
    
    if ~isempty(input.distance)
        ev = ev([ev.D]==input.distance); % narrow down to the right distance category
    end
    
    hold_state = input.axes.NextPlot;
    on_clean = onCleanup(@() set(input.axes, 'NextPlot', hold_state)); 
    
    data_failed = {};     
    data_passed = {};
    
    for ii = 1:length(input.pars)
        data_failed{ii} = [ev([ev.passed]==0).(input.pars(ii))];
        data_passed{ii} = [ev([ev.passed]==1).(input.pars(ii))];
    end
    
    data_failed{end+1} = input.marker_size/2;
    data_failed{end+1} = [ev([ev.passed]==0).detect_snr]; 
    data_passed{end+1} = input.marker_size;
    data_passed{end+1} = [ev([ev.passed]==1).detect_snr]; 
    
    if length(input.pars)==2 
        plot_func = @scatter; 
    elseif length(input.pars)>3
        plot_func = @scatter3; 
    else
        error('input "pars" must have 2-3 characters... got %s instead.', input.pars); 
    end
    
    h = plot_func(input.axes, data_failed{:}, '.'); 
    hold(input.axes, 'on'); 
    h = plot_func(input.axes, data_passed{:}, 'filled'); 

    input.axes.FontSize = input.font_size; 
    
    labels.R = 'stellar radius [FSU]'; 
    labels.r = 'occulter radius [FSU]'; 
    labels.b = 'impact parameter [FSU]'; 
    labels.v = 'transverse velocity [FSU/s]'; 
    
    xlabel(input.axes, labels.(input.pars(1))); 
    ylabel(input.axes, labels.(input.pars(2))); 
    
    colorbar(input.axes); 
    box(input.axes, 'on'); 
    
end








