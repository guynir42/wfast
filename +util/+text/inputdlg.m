function [R] = inputdlg(str,value)
% Usage: inputdlg(str, value=5)
% Displays an input dialog with question text "str" and default value 5.
% User may hit return after entering info.
% 
% taken shamelessly from https://www.mathworks.com/matlabcentral/newsreader/view_thread/295157
% and also added a cancel button

if nargin==0, help('util.text.inputdlg'); return; end

if nargin<2 || isempty(value)
    value = '5';
end

if ~ischar(value)
    value = num2str(value);
end

R = []; % In case the user closes the GUI.
S.fh = figure('units','pixels',...
              'position',[50 500 200 100],...
              'menubar','none',...
              'numbertitle','off',...
              'name',str,...
              'resize','off');
S.ed = uicontrol('style','edit',...
                 'units','pix',...
                'position',[10 60 180 30],...
                'string',value);
S.pb = uicontrol('style','pushbutton',...
                 'units','pix',...
                'position',[10 20 90 30],...
                'string','Enter',...
                'callback',{@pb_call});
S.cn = uicontrol('Style', 'pushbutton',...
                'units', 'pix',...
                'position',[100 20 90 30],...
                'string', 'Cancel',...
                'callback',{@cn_call});
set(S.ed,'call',@ed_call)
uicontrol(S.ed) % Make the editbox active.
uiwait(S.fh) % Prevent all other processes from starting until closed.

    function [] = pb_call(varargin)
%         R = str2double(get(S.ed,'string'));
        R = get(S.ed,'String');
        close(S.fh); % Closes the GUI, allows the new R to be returned.
    end

    function [] = ed_call(varargin)
        uicontrol(S.pb)
        drawnow
%         R = str2double(get(S.ed,'string'));
        if isvalid(S.ed)
            R = get(S.ed, 'String'); 
        else
            R = [];
        end
        close(gcbf)
    end
    
    function [] = cn_call(varargin)
        R = [];
        close(gcbf);
    end

end