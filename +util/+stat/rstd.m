function Rstd = rstd(Mat,Dim)
% Usage: Rstd = rstd(Mat,Dim)
% Calculate a robust estimator for the standard deviation of Mat along the  
% dimention Dim. 
% Will get the 25% and 75% percentile and calculate the average between them, 
% and multiply by 1.4826 to make it match the regular STD for a gaussian 
% noise source. 
%
% Taken shamelessly from Eran's code of similar name: Util.stat.rstd

    if nargin==0, help('util.stat.rstd'); return; end

    if nargin<2 || isempty(Dim)
        Dim = 1;
    end

    Factor = 1.4826;  % = 1./norminv(0.75,0,1)

    % think about making alternative notation that can handle Dim being a vector, or else just write a new rstd2 function for images... 
    ValLow  = prctile(Mat,25,Dim);
    ValHigh = prctile(Mat,75,Dim);
    Rstd    = (ValHigh - ValLow).*0.5.*Factor;
    
end