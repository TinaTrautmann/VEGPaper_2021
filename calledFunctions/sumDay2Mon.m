function [Out,NaNCnt] = sumDay2Mon(In,startDate,endDate)
% Calculates monthly nanmean
%
% Usages:
%   [Out,NaNCnt] = sumDay2Mon(In,startDate,endDate)
%
% Requires:
%   - In: daily time series of a variable (pix,tix)
%   - startDate: start date of time series 'YYYY-MM-DD'
%   - endDate: end date of time series 'YYYY-MM-DD'
%   - Out: time series with monthly sums (pix,tix) excluding NaNs
%   - NaNCnt: count of NaNs in each month (pix,tix)
%
% Purposes:
%
% Conventions:
%
% Created by:
%   - Tina Trautmann (ttraut)
%
% References:
%
%
% Versions:
%   - 2.0 on 27.09.2021


[xDay, xdays, xmonths,xyears, Y, M] = createDateVector(startDate,endDate,'d');

startYear   = min(Y);
endYear     = max(Y);


pix     = size(In,1);

Out     = NaN(pix,xmonths);
NaNCnt  = NaN(pix,xmonths);

cnt = 1;
for year=startYear:endYear
    
    for month=1:12
        % only average if there are valids
        valid = Y==year & M==month;
        if sum(valid(:)) ~= 0
            Out(:,cnt)      = nansum(In(:,valid),2);
            NaNCnt(:,cnt)   = sum(isnan(In(:,valid)),2);
            cnt = cnt+1;
        end
    end
    
end


end

