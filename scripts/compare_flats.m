%% loads a couple of calibration objects and show the different flats

% old calibration, before installing skirt adjustment
cal1 = img.Calibration;
cal1.camera_name = 'Balor'; 
cal1.loadByDate('2020-09-01');  

% after adding skirt adjustment
cal2 = img.Calibration;
cal2.camera_name = 'Balor'; 
cal2.loadByDate('2020-10-24');  

%%


f1 = util.plot.FigHandler('flat comparison'); 
f1.clear;

ax1 = axes('Parent', f1.fig, 'Position', [0.03 0.03 0.45 0.94]); 

util.plot.show(100*(cal1.flat_field-1), 'fancy', 'off', 'autodyn', 'on'); 

ax1.FontSize = 18;

ax2 = axes('Parent', f1.fig, 'Position', [0.5 0.03 0.45 0.94]); 

util.plot.show(100*(cal2.flat_field-1), 'fancy', 'off', 'autodyn', 'on'); 

ax2.CLim = ax1.CLim;

ax2.FontSize = ax1.FontSize; 













