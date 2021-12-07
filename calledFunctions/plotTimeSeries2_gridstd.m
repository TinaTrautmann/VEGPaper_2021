function [sname] = plotTimeSeries2_gridstd(name, tSteps, time, legendtext, In1, In2)
% Plots time series and MSC of up to 2 variables their std, areaweighted mean and median, returns save name
% name          - for title
% tSteps        - time steps of input data ('d' or 'm')
% time          - startDate & endDate
% legendtext    - cell structure with strings for legend (In1, In2)
% In1  - data 1 (red)
% In2  - data 2 (blue)




name        = strrep(name,'_',' ');
startDate   = time{1};
endDate     = time{2};
if tSteps == 'd'
    xTime      = createDateVector(startDate,endDate ,'m');
    [~,vecMonth]    = datevec(xTime);
    try
        dataM1           = aggDay2Mon(In1,startDate,endDate);
    end
    try
        dataM2           = aggDay2Mon(In2,startDate,endDate);
    end


elseif tSteps=='m'
    xTime      = createDateVector(startDate,endDate ,'m');
    [~,vecMonth]    = datevec(xTime);
    try
        dataM1           = In1;
    end
    try
        dataM2           = In2;        
    end

else
    error('no valid timestep given!')
end


% 
xSeason = [1:1:12];

% calc MSC + calc nanmean, nanmedian and std
if exist('dataM1', 'var') && ~isempty(dataM1)
    mean1   = mean(dataM1,1, 'omitnan');
    [~, meanSeason1] = calcMSC(dataM1,vecMonth);
    
    if size(dataM1,1)>1
        dataSeason1 = NaN(size(dataM1,1),12);
        for npix=1:size(dataM1,1)
            dataSeason1 (npix,:) = calcMSC(dataM1(npix,:),vecMonth);
        end
        
        median1 = median(dataM1,1, 'omitnan');
        medianSeason1   = median(dataSeason1,1, 'omitnan');
        
        std1        = std(dataM1,[],1,'omitnan');
        stdSeason1  = std(dataSeason1,[],1,'omitnan');
    end
end

if exist('dataM2', 'var') && ~isempty(dataM2)
    mean2   = mean(dataM2,1, 'omitnan');
    
    [~, meanSeason2] = calcMSC(dataM2,vecMonth);
    
    
    if size(dataM2,1)>1
        dataSeason2 = NaN(size(dataM2,1),12);
        for npix=1:size(dataM2,1)
            dataSeason2 (npix,:) = calcMSC(dataM2(npix,:),vecMonth);
        end
        
        median2 = median(dataM2,1, 'omitnan');
        medianSeason2   = median(dataSeason2,1, 'omitnan');
        
        std2        = std(dataM2,[],1,'omitnan');
        stdSeason2  = std(dataSeason2,[],1,'omitnan');
    end
end


% % correlation with first input data (daily values)
% try
%     cor_In3 = round(corr(In1(:),In3(:), 'rows', 'complete'),2);
% end

% include correlation in legendtext if In1 and In3 exists
if ~isempty(In1)
    if exist('std1', 'var')
        leg1 = {[legendtext{1} ' std'], [legendtext{1} ' mean'], [legendtext{1} ' median']};
    else
        leg1 = {[legendtext{1} ' mean']};
    end
else
    leg1 = [];
end

if  ~isempty(In2)
    if exist('std2', 'var')
        leg2 = {[legendtext{2} ' std'], [legendtext{2} ' mean'], [legendtext{2} ' median']};
    else
        leg2 = {[legendtext{2} ' mean']};
    end
else
    leg2 = [];
end

leg = [leg1, leg2];
leg = strrep(leg, '_','-');

%% Plot
figure('Color', [1 1 1]);
set(gcf, 'Position', [5 5 20 7])

%% Time Series
sp1 = subplot(1,2,1);
if ~isempty(In1)
    if exist('std1', 'var')
        boundedline(datenum(xTime),mean1,std1, 'alpha', 'nan', 'gap', 'cmap', rgb('DarkRed'), 'transparency', 0.05)
        hold on,
        plot(datenum(xTime),median1, ':', 'color', rgb('DarkRed'))
    else
        plot(datenum(xTime),mean1, '-', 'color', rgb('DarkRed'))
    end
end
hold on
if~isempty(In2)
    if exist('std2', 'var')
        boundedline(datenum(xTime),mean2,std2, 'alpha', 'nan', 'gap', 'cmap', rgb('Blue'), 'transparency', 0.05)        
        hold on,
        plot(datenum(xTime),median2, ':', 'color', rgb('Blue'))
    else
        plot(datenum(xTime),mean2, '-', 'color', rgb('Blue'))
    end
end

hold on, plot([datenum(xTime(1)) datenum(xTime(end))],[0 0], '-', 'color', [.5 .5 .5]) 
% set axis
box on, grid on
title(['time series (' tSteps ')'], 'fontsize', 9, 'fontweight', 'b');

years    = unique(year(xTime));
dn_array = NaN(length(years)+1,1);
for yy=1:length(years)
    dn_array(yy,1) = datenum(['01-Jan-' num2str(years(yy))]);
end
dn_array(end,1) = datenum(['01-Jan-' num2str(years(yy)+1)]);
dn_array(1:2:end) = [];
set(gca,'xtick',dn_array)
datetick('x','yyyy', 'keepticks')

%% MSC
sp2 = subplot(1,2,2);
if ~isempty(In1)
    if exist('stdSeason1', 'var')
        boundedline(xSeason,meanSeason1,stdSeason1, 'alpha', 'nan', 'gap', 'cmap', rgb('DarkRed'), 'transparency', 0.05)
        hold on,
        plot(xSeason,medianSeason1, ':', 'color', rgb('DarkRed'))
    else
        plot(xSeason,meanSeason1, '-', 'color', rgb('DarkRed'))
    end
end
hold on
if ~isempty(In2)
    if exist('stdSeason2', 'var')
        boundedline(xSeason,meanSeason2,stdSeason2, 'alpha', 'nan', 'gap', 'cmap', rgb('Blue'), 'transparency', 0.05)
        hold on,
        plot(xSeason,medianSeason2, ':', 'color', rgb('Blue'))
    else
        plot(xSeason,meanSeason2, '-', 'color', rgb('Blue'))
    end
end

hold on, plot([0.5 12.5],[0 0], '-', 'color', [.5 .5 .5]) 

% set axis
set(gca, 'Xlim', [0.5 12.5], 'XTick', [1:1:12])
box on, grid on
title('mean seasonal' , 'fontsize', 9, 'fontweight', 'b');

legend(leg, 'Position', [.44 .025 .12 .07], 'box', 'off', 'Orientation', 'Horizontal');

% add overall title
set(sp1,'fontsize',8,'Position', [.04 .22 .58 .64])
set(sp2,'fontsize',8,'Position', [.68 .22 .29 .64])

annotation('textbox', 'String', name, 'fontsize', 10, 'fontweight', 'bold',...
    'edgecolor', 'none', 'Position', [0 0.94 1 0.07], 'HorizontalAlignment', 'center' );


% save name
name        = strrep(name, ',', '_');
name        = strrep(name, '&', '_');
name        = regexprep(name,'\W','');
sname       = name;

end