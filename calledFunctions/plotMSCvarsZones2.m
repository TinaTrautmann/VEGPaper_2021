%% Plot MSC of obs & 2 simulations for 3 variables globally & for different zones
function [sname, sp] = plotMSCvarsZones2(rowNames, varNames, dataIn, dataObs, dataCR,CRnames, zoneNames, zoneIdx, pix_a, unitNames, col)
% plots the global and zonal MSC of 1-2 experiments against observations
% for different variables
%
% returns the save name of the figure, the figure handle and a
% table-structure with the correlation of the 1-2 experiments with
% observations for all zones and variables
%
% Input:
%     rowNames      = experiment names
%     colNames      = variable names to plot
%     dataIn        = structure of variables, fields = variables with MSC for each grid; size(npix,12)
%     dataObs       = structure with observational data to plot, if it is
%     an input, the correlation with obs is calculated instead of the
%     contribution ratio
%     zoneNames     = name of the zones defined, order should reflect the
%                       value in zoneIdx (i.e. first name refers to
%                       zoneIdx = 1
%     zoneIdx       = classifies the grids to different zones, continous
%                       numbers starting from 1; size(npix,1)
%     pix_a         = area of each grid; size(npix,1)
%     unitnames     = cellstring with units of each variable, can be empty
%     col           = line color for each variable that's plotted
%
% Output:
%     sname     = save name (= title of the figure)
%     T_corr    = correlation table for each variable, combined as a
%                   structure
%     sp        = figure handle

%% checks?

%%
if isempty(col)
    col         = othercolor('Paired10',numel(varNames))
end

expNames  = strrep(rowNames,'_','-');
zoneNames = strrep(zoneNames,'_','-');

expNamesTmp     = expNames;
try
idx             = strcmp(expNamesTmp,'E-B-bL-RD4')
expNamesTmp{idx==1} = 'VEG';
end

zNames    = ['Global', zoneNames];
zNames2   = strrep(zNames,'-','');
zNames2   = strrep(zNames2,' ','');

nPix    = length(pix_a);
nrows   = numel(rowNames);
ncols   = numel(zoneNames)+1;
xSeason = 1:1:12;

sp  = figure; set(gcf, 'Position', [5 2 ncols*5 nrows*5]);
if nrows==1
    set(gcf, 'Position', [5 1 ncols*5 nrows*5+2]);
    ha  = tight_subplot(nrows,ncols,[.05 .02],[.15 .15],[.05 .02]);
else
    set(gcf, 'Position', [5 1 ncols*5 nrows*5]);
    ha  = tight_subplot(nrows,ncols,[.05 .02],[.08 .08],[.05 .02]);
end

% loop over rows (=experiments), then cols (=global + zones)
cnt = 1;
for rr=1:nrows
    expN  = rowNames{rr};
    
    % loop of zones
    for cc=1:ncols
        if cc==1
            idxZ = 1:1:nPix;
        else
            idxZ = find(zoneIdx==cc-1);
        end
        axes(ha(cnt));
        
        if ~isempty(dataObs)
            obsN = fieldnames(dataObs);
            obsN = obsN{1};
            zData.(obsN) = nanmeanArea(dataObs.(obsN)(idxZ,:),pix_a(idxZ));
            plot(xSeason,zData.(obsN), 'LineStyle', ':', 'LineWidth', 1.5, 'Color', rgb('DarkRed')), hold on
        end
        
        %loop over variables
        for vN=1:numel(varNames)
            varN = varNames{vN};
            zData.(varN) = nanmeanArea(dataIn.(expN).(varN)(idxZ,:),pix_a(idxZ));
            plot(xSeason,zData.(varN), 'LineStyle', '-', 'Color', col(vN,:)), hold on
        end
        plot([0.5 12.5],[0 0], '-', 'color', [.5 .5 .5])
        if rr==1; title(zNames{cc}, 'Fontsize',10); end
        if cc==1; yl= ylabel(expNamesTmp{rr}, 'Fontweight', 'b'); yl.Position(1) = -1.5; end
        ha(cnt).XAxis.MinorTickValues =  [1:1:12];
        
        
        % correlation if there are observational data
        if ~isempty(dataObs)
            for vN=1:numel(varNames)
                varN = varNames{vN};
                T_corr.(varN) = corr(zData.(obsN)(:),zData.(varN)(:),'Rows', 'complete');
%                 if cc==5
%                     t1 = text(0.97,0+0.05*vN, ['corr ' varN ' = ' num2str(round(T_corr.(varN),2)) ],'Units','Normalized', 'Fontsize', 6, 'Color', col(vN,:), 'HorizontalAlignment', 'right');
%                 else
%                     t1 = text(0.03,0+0.05*vN, ['corr ' varN ' = ' num2str(round(T_corr.(varN),2)) ],'Units','Normalized', 'Fontsize', 6, 'Color', col(vN,:), 'HorizontalAlignment', 'left');
%                 end
                if cc==5
                    t1 = text(0.03,1-0.05*vN, ['corr ' varN  ' = ' num2str(round(T_corr.(varN),2)) ],'Units','Normalized', 'Fontsize', 6, 'Color', col(vN,:), 'HorizontalAlignment', 'left');
                else
                    t1 = text(0.97,1-0.05*vN, ['corr ' varN ' = ' num2str(round(T_corr.(varN),2)) ],'Units','Normalized', 'Fontsize', 6, 'Color', col(vN,:), 'HorizontalAlignment', 'right');
                end
            end
        elseif ~isempty(dataCR) % calculate the CR for each zone, if there are only 3 input
            % variables
            iNames = CRnames;
            for vN=1:numel(iNames)
                iName   = iNames{vN};
                tmpData = nanmeanArea(dataCR.(expN).(iName)(idxZ,:),pix_a(idxZ));
                iName   = strrep(iName, 'I', 'I-');
                if cc==5
                    t1 = text(0.97,0+0.07*vN, [iName ' = ' num2str(round(tmpData,2)) ],'Units','Normalized', 'Fontsize', 8, 'Color', col(vN+1,:), 'HorizontalAlignment', 'right');
                else
                    t1 = text(0.03,0+0.07*vN, [iName ' = ' num2str(round(tmpData,2)) ],'Units','Normalized', 'Fontsize', 8, 'Color', col(vN+1,:), 'HorizontalAlignment', 'left');
                end
            end
        end
        
        cnt = cnt+1;
        
    end
end
% set axis
set(ha(:), 'Xlim', [0.5 12.5], 'XTick', [3:3:12], 'GridAlpha', 0.25, 'XMinorGrid', 'on', 'MinorGridAlpha', 0.05, 'MinorGridLineStyle', '-',...
    'YGrid', 'on', 'Fontsize', 8)
set(ha(1:end-ncols), 'XTickLabel', [])

for tmp=1:length(ha),ha(tmp).YLabel.FontSize = 10;end

% legend & overall title
varNames_org = varNames;
varNames = strrep(varNames,'I','I-');
varNames = strrep(varNames,'wGW','wDeep');
varNames = strrep(varNames,'wSurf','wSlow');

if ~isempty(dataObs)
    obsN = fieldnames(dataObs);
    varNames = [obsN, varNames_org'];
end
l  = legend(strrep(varNames,'_','-'), 'Position', [.44 -.01 .12 .07], 'box', 'off', 'Orientation', 'Horizontal','Fontsize', 8);

if nrows==1
    l.Position =  [.44 0 .12 .07];
else                    
    l.Position =  [.44 -.02 .11 .07];
end
an = annotation('textbox',[0 .52 1 .5],'String', ['TWS Composition - Mean Seasonal Cycle' ] ,'LineStyle','none', 'Fontsize', 10, 'Fontweight', 'b','fitboxtotext', 'on', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center');

% savename
sname = ['TWS Composition -  Mean Seasonal Cycle - ' strjoin(strrep(varNames,'_','-'), ' vs ') ];

end