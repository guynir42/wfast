% This script has some examples for running the camera and acquisition class
% from the command line. Make sure to input your own set of parameters! 

%% load the Acquisition object along with the camera. 

if ~exist('a', 'var') || isempty(a) || ~isa(a, 'img.Acquisition') % verify this object is not already loaded! 
    a = img.Acquisition; a.cal.load; a.chooseSource('camera'); % create an object, load the calibration, choose the camera as source
end

%% Use the default focus position

a.cam.setupDefaultFocusPosition; 

%% Use the autofocus module of Acquisition

a.runFocus;

%% use the camera to take full frame images

a.cam.run('use_save', 1, 'expT', 3, 'frame_rate', NaN, 'num_batches', 10, 'batch_size', 1); % use this to take 10 images at 3 second exposure, as fast as it can, one image per file. 

a.cam.run('use_save', 1, 'expT', 0.025, 'frame_rate', 30, 'num_batches', 10, 'batch_size', 100); % use this to take 10 batches of 100 images, at 0.025 second, at constant frame rate. 

%% use the pipeline of Acquisition class to produce cutouts / stacks

a.run('use_save', 1, 'expT', 0.025, 'frame_rate', 25, 'num_batches', 10, 'batch_size', 100, 'num_stars', 750, 'cut_size', 21, 'num_backgrounds', 50, 'cut_size_bg', 30, ...
    'objname', 'KBO_survey', 'RA', '11 22 33.44', 'DE', '+30 40 50'); % the full list of inputs is in Acquisition/makeInputVars

% NOTE: the inputs used when calling run(...) are not saved in the class. 
%       To change the object default use, e.g., a.expT = 3; 

%% deflate your data:

a.deflator.run('src', 'C:\data_temp\2019-06-01', 'out', 'D:\Dropbox\data\WFAST\2019-06-01'); % make sure to change to the relevant directory. Out dir will be created. 

% NOTE: make sure to delete the data from C:\ after it has been copied/deflated. 

