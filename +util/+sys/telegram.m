function ret = telegram(varargin)
% Usage: ret = telegram(varargin)
% Send a message to a Telegram bot. 
%
% Use telegram(varargin) in the same way as sprintf(varargin)
% Example: telegram('Hello, World!');
%          telegram('%d + %d = %d',1,2,1+2);
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
% Getting the tokens: https://solvit.io/0f9c61a

    if nargin==0; help('util.sys.telegram'); return; end

    % default token and chat_id
    token = '978929699:AAE9m1EITSSWdfJMVgWHAuI0gJnlEku1aZ0';
    chat_id = '1121382138';

    str = sprintf(varargin{:});

    % print to MATLAB command window
    fprintf(str);

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
