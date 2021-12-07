%% Plot MSC of obs & 2 simulations for 3 variables globally & for different zones
function [sname, sp] = plotTSvarsZones2(figNames, varNames, dataIn, dataObs, dataCR, CRnames, zoneNames, zoneIdx, pix_a, xTime, col, spth)
% plots the global and zonal time series of n experiments
% for different variables
%
% returns the save name of the figure, the figure handle and a
% table-structure with the correlation of the 1-2 experiments with
% observations for all zones and variables
%
% Input:
%     figNames      = experiment names
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

expNames  = strrep(figNames,'_','-');
zoneNames = strrep(zoneNames,'_','-');


zNames    = ['Global', zoneNames];

nPix    = length(pix_a);
nrows   = numel(zoneNames)+1;
ncols   = 1;

% xTime      = createDateVector(startDate,endDate ,'m');

%% loop over experiments
for ff=1:numel(figNames)
    expName = figNames{ff};
    expName2 = strrep(expName,'E_B_bL_RD4', 'VEG');
    sp  = figure; set(gcf, 'Position', [5 2 15 nrows*4]);
    ha  = tight_subplot(nrows,ncols,[.03 .02],[.05 .05],[.08 .02]);
    % loop over rows (=experiments), then cols (=global + zones)
    for rr=1:nrows
        expN  = zNames{rr};
        
        if rr==1
            idxZ = 1:1:nPix;
        else
            idxZ = find(zoneIdx==rr-1);
        end
        axes(ha(rr));
        
        if ~isempty(dataObs)
            obsN = fieldnames(dataObs);
            obsN = obsN{1};
            zData.(obsN) = nanmeanArea(dataObs.(obsN)(idxZ,:),pix_a(idxZ));
            plot(xTime,zData.(obsN), 'LineStyle', ':', 'LineWidth', 1.5, 'Color', rgb('DarkRed')), hold on
        end
        
        %loop over variables
        for vN=1:numel(varNames)
            varN = varNames{vN};
            zData.(varN) = nanmeanArea(dataIn.(expName).(varN)(idxZ,:),pix_a(idxZ));
            
            plot(xTime,zData.(varN), 'LineStyle', '-', 'Color', col(vN,:)), hold on
        end
        plot([xTime(1) xTime(end)],[0 0], '-', 'color', [.5 .5 .5])
        yl= ylabel(zNames{rr}, 'Fontweight', 'b'); yl.Position(1) = -200;
        
        % correlation if there are observational data
        if ~isempty(dataObs)
            for vN=1:numel(varNames)
                varN = varNames{vN};
                T_corr.(varN) = corr(zData.(obsN)(:),zData.(varN)(:),'Rows', 'complete');
                t1 = text(0.01,1-0.07*vN, ['corr ' strrep(varN,'_','-') ' = ' num2str(round(T_corr.(varN),2)) ],'Units','Normalized', 'Fontsize', 6, 'Color', col(vN,:), 'HorizontalAlignment', 'left');
            end
        elseif ~isempty(dataCR) % calculate the CR for each zone, if there are only 3 input
            % variables
            iNames = CRnames;
            for vN=1:numel(iNames)
                iName   = iNames{vN};
                tmpData = nanmeanArea(dataCR.(expName).(iName)(idxZ,:),pix_a(idxZ));
                iName   = strrep(iName, 'I', 'I-');
                if rr==1 || rr==3 || rr==5 %left up
                    t1 = text(0.03,0.6+0.07*vN, [iName ' = ' num2str(round(tmpData,2)) ],'Units','Normalized', 'Fontsize', 6, 'Color', col(1+vN,:), 'HorizontalAlignment', 'left');
                else %left bottom
                    t1 = text(0.03,0+0.07*vN, [iName ' = ' num2str(round(tmpData,2)) ],'Units','Normalized', 'Fontsize', 6, 'Color', col(1+vN,:), 'HorizontalAlignment', 'left');
                end
            end
            
        end
        % set axis
        set(ha(:), 'GridAlpha', 0.25, 'XMinorGrid', 'on', 'MinorGridAlpha', 0.05, 'MinorGridLineStyle', '-',...
            'YGrid', 'on', 'Fontsize', 7)
        
    end
    for tmp=1:length(ha),ha(tmp).YLabel.FontSize = 10;end
    
    % legend & overall title
    if ~isempty(dataObs)
        obsN = fieldnames(dataObs);
        legNames = [obsN, varNames];
    else
        legNames = varNames;
    end
    
    legNames = strrep(legNames,'wGW','wDeep');
    legNames = strrep(legNames,'wSurf','wSlow');
    
    l  = legend(strrep(legNames,'_','-'), 'Position', [.44 -.02 .12 .07], 'box', 'off', 'Orientation', 'Horizontal','Fontsize', 8);
    an = annotation('textbox',[0 .5 1 .5],'String', ['TWS Composition - IAV | ' strrep(expName2,'_','-') ' | ' strjoin(strrep(varNames,'_','-'), ' vs ') ] ,'LineStyle','none', 'Fontsize', 10, 'Fontweight', 'b','fitboxtotext', 'on', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center');
    
    % savename
    sname = ['TWS Composition -  IAV - ' strrep(expName2,'_','-') ' - ' strjoin(strrep(varNames,'_','-'), ' vs ') ];
    % print
    print(gcf,[spth sname '.png'],'-dpng','-r300');
    
end

end