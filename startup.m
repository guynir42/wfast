addpath(getenv('WFAST'));
addpath(getenv('MAAT'));

if ~strcmp(version('-release'), '2017a')
    set(groot,'defaultFigureCreateFcn',@(fig,~)addToolbarExplorationButtons(fig))
    set(groot,'defaultAxesCreateFcn',@(ax,~)set(ax.Toolbar,'Visible','off'))
end
