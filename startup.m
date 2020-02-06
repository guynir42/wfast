addpath(getenv('WFAST'));
addpath(getenv('MAAT'));

v = version('-release');

if ~strcmp(v, '2017a')
    set(groot,'defaultFigureCreateFcn',@(fig,~)addToolbarExplorationButtons(fig))
    set(groot,'defaultAxesCreateFcn',@(ax,~)set(ax.Toolbar,'Visible','off'))
end

