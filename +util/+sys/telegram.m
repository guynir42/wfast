function ret = telegram(token, chat_id, varargin)
% Usage: ret = telegram(token, chat_id, varargin)
% Send a message using a Telegram bot.
% Inputs: -token is the bot identifier you got when the bot was created.
%         -chat_id is the ID of the user you want to send to.
%
% Use telegram(token, chat_id, varargin) in the same way as sprintf(varargin)
% Example: telegram(token, chat_id, 'Hello, World!');
%          telegram(token, chat_id, '%d + %d = %d',1,2,1+2);
% 
% Define token and chat_id before use, 
% which are the authorization token of the target Telegram bot 
% and the identifier or username of the target chat
%
% Please refer the following post 
% "Creating a Telegram bot for personal notifications"
% https://www.forsomedefinition.com/automation/creating-telegram-bot-notifications/
% 
% Seongsik Park
% seongsikpark@postech.ac.kr
%
% Getting the token and chat_id: https://solvit.io/0f9c61a

    if nargin==0; help('util.sys.telegram'); return; end

    if nargin<3 || isempty(token) || isempty(chat_id)
        error('Must supply a bot token and a chat_id for this to work!');
    end
    
    str = sprintf(varargin{:});

    % print to MATLAB command window
%     fprintf(str);

    % convert MATLAB string to url query string
    sendstr = urlencode(str);
    sendstr = ['https://api.telegram.org/bot',token,...
               '/sendMessage?chat_id=',chat_id,...
               '&text=',sendstr];

    % send a message   
    ret = webread(sendstr); 
    assert(ret.ok);

    % append human readable datetime to results [Set TimeZone value to desired time zone]
    ret.result.datetime=datetime(ret.result.date,'ConvertFrom','posixtime','TimeZone','Asia/Seoul');

end
