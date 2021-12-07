function [ idx_1, idx_2, time3, xData3] = calcConsisTime( time1, time2, timesteps )
% [ idx_1, idx_2, time3, xData3] = CalcConsisTime( time1, time2, timesteps )
% returns index of consisting time period of 2 input time periods, in resolution
% defined by timesteps (monthly or daily)

if verLessThan('matlab','8.4') % 8.4 is R2014b -> needed for datetime function
    msg = 'Matlab Version R2014b or newer needed to perform calcConsisTime!'
    error(msg)
end

try
    % Time1
    startD1      = datetime(time1(1),'InputFormat','MM-dd-yyyy');
    endD1        = datetime(time1(2),'InputFormat','MM-dd-yyyy');
    
    % Time2
    startD2      = datetime(time2(1),'InputFormat','MM-dd-yyyy');
    endD2        = datetime(time2(2),'InputFormat','MM-dd-yyyy');
catch
    % Time1
    startD1      = datetime(time1(1),'InputFormat','yyyy-MM-dd');
    endD1        = datetime(time1(2),'InputFormat','yyyy-MM-dd');
    
    % Time2
    startD2      = datetime(time2(1),'InputFormat','yyyy-MM-dd');
    endD2        = datetime(time2(2),'InputFormat','yyyy-MM-dd');
end

% monthly or daily timesteps required?
switch timesteps
    case 'daily'
        % Time1
        xDay1        = [startD1:endD1];
        [Y1,M1,D1]   = datevec(xDay1);
        % Time2
        xDay2        = [startD2:endD2];
        [Y2,M2,D2]   = datevec(xDay2);
        
        % find idx (cols) that coincidence
        idx_1   = find(ismember(xDay1,xDay2)); % cols of data1 that are in data2
        idx_2   = find(ismember(xDay2,xDay1)); % cols of data2 that are in data1

%         idx_1   = find(ismember(Y1,Y2) & ismember(M1,M2) & ismember(D1,D2)); % cols of data1 that are in data2
%         idx_2   = find(ismember(Y2,Y1) & ismember(M2,M1) & ismember(D2,D1)); % cols of data2 that are in data1
       
        if isempty(idx_1) == 1 || isempty(idx_2) ==1
            error('Timeperiods do not overlap!')
        end
        
        % consistent time period
        xData3      = xDay2(idx_2);
        [Y3,M3,D3]  = datevec(xData3);


%%----------------------------------------------------------------------------
    case 'monthly'
        % Time1
        months1      = calmonths(between(startD1, endD1,'months'))+1; %+1 for column of start month
        xMonth1      = [startD1,startD1+calmonths(1:months1-1)];
        
        % Time2
        months2      = calmonths(between(startD2, endD2,'months'))+1; %+1 for column of start month
        xMonth2      = [startD2,startD2+calmonths(1:months2-1)];
        
        % find idx (cols) that coincidence
        idx_1   = find(ismember(xMonth1,xMonth2)); % cols of data1 that are in data2
        idx_2   = find(ismember(xMonth2,xMonth1)); % cols of data2 that are in data1
        
        if isempty(idx_1) == 1 || isempty(idx_2) ==1
            error('Timeperiods do not overlap!')
        end
        
        % consistent time period
        xData3      = xMonth2(idx_2);
        [Y3,M3,D3]  = datevec(xData3);
        D3(end)     = eomday(Y3(end),M3(end));
        xData3(end) = datestr(datenum(Y3(end), M3(end), D3(end)), 'dd-mmm-yyyy');
%% ----------------------------------------------------------------------------
    otherwise
        error(['Define timesteps!' char(10) ...
            'It has to be either' char(39) 'daily' char(39) ' or ' char(39) 'monthly' char(39) '!'])
        
end

% output
startDate3  = datestr(datenum(Y3(1), M3(1), D3(1)), 'mm-dd-yyyy');
endDate3    = datestr(datenum(Y3(end), M3(end), D3(end)), 'mm-dd-yyyy');
time3       = {startDate3, endDate3};


end

