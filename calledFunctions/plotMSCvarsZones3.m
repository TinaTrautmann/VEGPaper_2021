%% Plot MSC of obs & 3 simulations for 3 variables globally & for different zones
function [sname, T_corr, T_AIC, T_BIC, sp] = plotMSCvarsZones3(dataNames, dataObs, data1, data2, data3, n_params, zoneNames, zoneIdx, pix_a, unitNames, col)
% plots the global and zonal MSC of 1-2 experiments against observations
% for different variables
%
% returns the save name of the figure, the figure handle and a
% table-structure with the correlation, AIC & BIC of the 1-3 experiments with
% observations for all zones and variables
%
% Input:
%     dataNames     = experiment names, in the order dataObs, data1, data2
%     dataObs       = structure of observations, fields = variables with MSC for each grid; size(npix,12)
%     data1         = structure of experiment1, fields = variables with MSC for each grid; size(npix,12)
%     data2         = structure of experiment2, fields = variables with MSC for each grid; size(npix,12)
%     data3         = structure of experiment3, fields = variables with MSC for each grid; size(npix,12)
%     zoneNames     = name of the zones defined, order should reflect the
%                       value in zoneIdx (i.e. first name refers to
%                       zoneIdx = 1
%     zoneIdx       = classifies the grids to different zones, continous
%                       numbers starting from 1; size(npix,1)
%     pix_a         = area of each grid; size(npix,1)
%     unitnames     = cellstring with units of each variable, can be empty
%     col           = line color for plotting obs, experiment1 and
%     experiment 2; default is red, blue, green
%
% Output:
%     sname     = save name (= title of the figure)
%     T_corr    = correlation table for each variable, combined as a
%                   structure
%     T_AIC    = correlation table for each variable, combined as a
%                   structure
%     T_BIC    = correlation table for each variable, combined as a
%                   structure
%     sp        = figure handle

% % ----EXAMPLE Input----
% dataNames = {'obs', 'exp1', 'exp2'};
% dataObs.A = randi(5,100,12);
% data1.A   = dataObs.A  ./2 ;
% data2.A   = dataObs.A  + 1.3;
% 
% dataObs.B = randi(5,100,12);
% data1.B   = (dataObs.A ./2 ) + 1;
% data2.B   =  data1.A + 1.2;
% 
% dataObs.C = randi(5,100,12);
% data1.C   = (dataObs.A ./2 ) + 1;
% data2.C   =  data1.A + 1.2;
% 
% zoneNames   = {'zone1','zone2','zone3','zone4','zone5'};
% zoneIdx     = randi(5,100,1);
% 
% pix_a       = ones(100,1);
% unitNames   = {'mm', 'mm/d', 'mm/d'};

%% checks?

%%
if isempty(col)
    col         = [rgb('DarkRed'); rgb('Blue'); rgb('Green'); rgb('DarkOrchid')];
end

dataNames = strrep(dataNames,'_','-');
zoneNames = strrep(zoneNames,'_','-');


varNames  = fieldnames(dataObs);
zNames    = ['Global', zoneNames];
zNames2   = strrep(zNames,'-','');
zNames2   = strrep(zNames2,' ','');

nrows   = numel(varNames);
ncols   = numel(zoneNames)+1;
xSeason = 1:1:12;

% correlations -> 3 tables, 1 per variable
tmpArray    = NaN(numel(dataNames)-1,numel(zNames));
varTypes    = cell(numel(zNames),1);
varTypes(:) = {'double'};
for vn=1:numel(varNames)
    T_corr.(varNames{vn}) = array2table(tmpArray,'VariableNames', zNames2, 'RowNames', dataNames(2:end));
    T_AIC.(varNames{vn}) = array2table(tmpArray,'VariableNames', zNames2, 'RowNames', dataNames(2:end));
    T_BIC.(varNames{vn}) = array2table(tmpArray,'VariableNames', zNames2, 'RowNames', dataNames(2:end));
end

sp  = figure; set(gcf, 'Position', [5 5 ncols*5 nrows*5]);
ha  = tight_subplot(nrows,ncols,[.05 .02],[.08 .08],[.05 .02]);
% loop over rows (=variables), then cols (=global + zones)
cnt = 1;
for rr=1:nrows
    varN  = varNames{rr};
    unitN = unitNames{rr};
    dObs = dataObs.(varN);
    d1   = data1.(varN);
    d2   = data2.(varN);
    d3   = data3.(varN);
    for cc=1:ncols
        if cc==1
            idxZ = 1:1:size(dObs,1);
        else
            idxZ = find(zoneIdx==cc-1);
        end
        axes(ha(cnt));
        plot(xSeason,nanmeanArea(dObs(idxZ,:),pix_a(idxZ)), 'Linewidth', 1.5, 'LineStyle', ':', 'Color', col(1,:)), hold on
        plot(xSeason,nanmeanArea(d1(idxZ,:),pix_a(idxZ)), 'LineStyle', '-', 'Color', col(2,:)), hold on
        plot(xSeason,nanmeanArea(d2(idxZ,:),pix_a(idxZ)), 'LineStyle', '-', 'Color', col(3,:)), hold on
        plot(xSeason,nanmeanArea(d3(idxZ,:),pix_a(idxZ)), 'LineStyle', '-', 'Color', col(4,:)), hold on
        plot([0.5 12.5],[0 0], '-', 'color', [.5 .5 .5])
        if rr==1; title(zNames{cc}); end
        if cc==1; yl= ylabel([varN ' [' unitN ']'], 'Fontweight', 'b'); yl.Position(1) = -1.5; end
        ha(cnt).XAxis.MinorTickValues =  [1:1:12];
        cnt = cnt+1;
        
        %calculate correlation
        tmpO = dObs(idxZ,:);
        tmp1 = d1(idxZ,:);
        tmp2 = d2(idxZ,:);
        tmp3 = d3(idxZ,:);
        T_corr.(varN){1,cc} = round(corr(tmpO(:),tmp1(:), 'rows', 'complete'),2);
        T_corr.(varN){2,cc} = round(corr(tmpO(:),tmp2(:), 'rows', 'complete'),2);
        T_corr.(varN){3,cc} = round(corr(tmpO(:),tmp3(:), 'rows', 'complete'),2);
        
        % calculate AIC & BIC
        aic_obs = tmpO(:);
        aic_mod = [tmp1(:), tmp2(:), tmp3(:)];
        [AIC, BIC] = calcAIC(aic_obs, aic_mod, n_params); 
        
        % delta AIC
        T_AIC.(varN){1,cc} = AIC(2,1);
        T_AIC.(varN){2,cc} = AIC(2,2);
        T_AIC.(varN){3,cc} = AIC(2,3);
        
        % delta BIC
        T_BIC.(varN){1,cc} = BIC(2,1);
        T_BIC.(varN){2,cc} = BIC(2,2);
        T_BIC.(varN){3,cc} = BIC(2,3);
        
        
        %add the correlation in the plot
        if strcmp(varN, 'wTWS')==1
            if cc==5
                t1 = text(0.97,0.275, ['r^2 = ' num2str(T_corr.(varN){1,cc}) ' | AIC = ' num2str(round(T_AIC.(varN){1,cc},2))],'Units','Normalized', 'Fontsize', 6, 'Color', col(2,:), 'HorizontalAlignment', 'right');
                t2 = text(0.97,0.175, ['r^2 = ' num2str(T_corr.(varN){2,cc}) ' | AIC = ' num2str(round(T_AIC.(varN){2,cc},2))],'Units','Normalized', 'Fontsize', 6, 'Color', col(3,:), 'HorizontalAlignment', 'right');
                t3 = text(0.97,0.075, ['r^2 = ' num2str(T_corr.(varN){3,cc}) ' | AIC = ' num2str(round(T_AIC.(varN){3,cc},2))],'Units','Normalized', 'Fontsize', 6, 'Color', col(4,:), 'HorizontalAlignment', 'right');
            else
                t1 = text(0.03,0.275, ['r^2 = ' num2str(T_corr.(varN){1,cc}) ' | AIC = ' num2str(round(T_AIC.(varN){1,cc},2))],'Units','Normalized', 'Fontsize', 6, 'Color', col(2,:), 'HorizontalAlignment', 'left');
                t2 = text(0.03,0.175, ['r^2 = ' num2str(T_corr.(varN){2,cc}) ' | AIC = ' num2str(round(T_AIC.(varN){2,cc},2))],'Units','Normalized', 'Fontsize', 6, 'Color', col(3,:), 'HorizontalAlignment', 'left');
                t3 = text(0.03,0.075, ['r^2 = ' num2str(T_corr.(varN){3,cc}) ' | AIC = ' num2str(round(T_AIC.(varN){3,cc},2))],'Units','Normalized', 'Fontsize', 6, 'Color', col(4,:), 'HorizontalAlignment', 'left');
            end
        elseif strcmp(varN, 'ET')==1
            if cc==2 || cc==3
                t1 = text(0.97,0.725, ['r^2 = ' num2str(T_corr.(varN){1,cc}) ' | AIC = ' num2str(round(T_AIC.(varN){1,cc},2))],'Units','Normalized', 'Fontsize', 6, 'Color', col(2,:), 'HorizontalAlignment', 'right');
                t2 = text(0.97,0.825, ['r^2 = ' num2str(T_corr.(varN){2,cc}) ' | AIC = ' num2str(round(T_AIC.(varN){2,cc},2))],'Units','Normalized', 'Fontsize', 6, 'Color', col(3,:), 'HorizontalAlignment', 'right');
                t3 = text(0.97,0.925, ['r^2 = ' num2str(T_corr.(varN){3,cc}) ' | AIC = ' num2str(round(T_AIC.(varN){3,cc},2))],'Units','Normalized', 'Fontsize', 6, 'Color', col(4,:), 'HorizontalAlignment', 'right');

            else
                t1 = text(0.97,0.275, ['r^2 = ' num2str(T_corr.(varN){1,cc}) ' | AIC = ' num2str(round(T_AIC.(varN){1,cc},2))],'Units','Normalized', 'Fontsize', 6, 'Color', col(2,:), 'HorizontalAlignment', 'right');
                t2 = text(0.97,0.175, ['r^2 = ' num2str(T_corr.(varN){2,cc}) ' | AIC = ' num2str(round(T_AIC.(varN){2,cc},2))],'Units','Normalized', 'Fontsize', 6, 'Color', col(3,:), 'HorizontalAlignment', 'right');
                t3 = text(0.97,0.075, ['r^2 = ' num2str(T_corr.(varN){3,cc}) ' | AIC = ' num2str(round(T_AIC.(varN){3,cc},2))],'Units','Normalized', 'Fontsize', 6, 'Color', col(4,:), 'HorizontalAlignment', 'right');
            end
        elseif strcmp(varN, 'Q')==1
            if cc==2 || cc==4
                t1 = text(0.97,0.725, ['r^2 = ' num2str(T_corr.(varN){1,cc}) ' | AIC = ' num2str(round(T_AIC.(varN){1,cc},2))],'Units','Normalized', 'Fontsize', 6, 'Color', col(2,:), 'HorizontalAlignment', 'right');
                t2 = text(0.97,0.825, ['r^2 = ' num2str(T_corr.(varN){2,cc}) ' | AIC = ' num2str(round(T_AIC.(varN){2,cc},2))],'Units','Normalized', 'Fontsize', 6, 'Color', col(3,:), 'HorizontalAlignment', 'right');
                t3 = text(0.97,0.925, ['r^2 = ' num2str(T_corr.(varN){3,cc}) ' | AIC = ' num2str(round(T_AIC.(varN){3,cc},2))],'Units','Normalized', 'Fontsize', 6, 'Color', col(4,:), 'HorizontalAlignment', 'right');
                
            else
                t1 = text(0.97,0.275, ['r^2 = ' num2str(T_corr.(varN){1,cc}) ' | AIC = ' num2str(round(T_AIC.(varN){1,cc},2))],'Units','Normalized', 'Fontsize', 6, 'Color', col(2,:), 'HorizontalAlignment', 'right');
                t2 = text(0.97,0.175, ['r^2 = ' num2str(T_corr.(varN){2,cc}) ' | AIC = ' num2str(round(T_AIC.(varN){2,cc},2))],'Units','Normalized', 'Fontsize', 6, 'Color', col(3,:), 'HorizontalAlignment', 'right');
                t3 = text(0.97,0.075, ['r^2 = ' num2str(T_corr.(varN){3,cc}) ' | AIC = ' num2str(round(T_AIC.(varN){3,cc},2))],'Units','Normalized', 'Fontsize', 6, 'Color', col(4,:), 'HorizontalAlignment', 'right');
            end
        else
            t1 = text(0.97,0.275, ['r^2 = ' num2str(T_corr.(varN){1,cc}) ' | AIC = ' num2str(round(T_AIC.(varN){1,cc},2))],'Units','Normalized', 'Fontsize', 6, 'Color', col(2,:), 'HorizontalAlignment', 'right');
            t2 = text(0.97,0.175, ['r^2 = ' num2str(T_corr.(varN){2,cc}) ' | AIC = ' num2str(round(T_AIC.(varN){2,cc},2))],'Units','Normalized', 'Fontsize', 6, 'Color', col(3,:), 'HorizontalAlignment', 'right');
            t3 = text(0.97,0.075, ['r^2 = ' num2str(T_corr.(varN){3,cc}) ' | AIC = ' num2str(round(T_AIC.(varN){3,cc},2))],'Units','Normalized', 'Fontsize', 6, 'Color', col(4,:), 'HorizontalAlignment', 'right');
        end
    end
end
% set axis
set(ha(:), 'Xlim', [0.5 12.5], 'XTick', [3:3:12], 'GridAlpha', 0.25, 'XMinorGrid', 'on', 'MinorGridAlpha', 0.05, 'MinorGridLineStyle', '-',...
    'YGrid', 'on', 'Fontsize', 7)
set(ha(1:end-ncols), 'XTickLabel', [])

for tmp=1:length(ha),ha(tmp).YLabel.FontSize = 10;end

% legend & overall title
tmpT = dataNames;

idx = strcmp(tmpT,'E-B-bL-RD4');
try
tmpT{idx==1} = 'VEG';
end
                    
l  = legend(tmpT, 'Position', [.44 -.02 .12 .07], 'box', 'off', 'Orientation', 'Horizontal','Fontsize', 7);
an = annotation('textbox',[0 .5 1 .5],'String', ['Comparison of the Mean Seasonal Cycle - ' strjoin(tmpT, ' vs ') ] ,'LineStyle','none', 'Fontsize', 10, 'Fontweight', 'b','fitboxtotext', 'on', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center');

% savename
sname = ['Comparison of the Mean Seasonal Cycle - ' strjoin(tmpT, ' vs ') ];

end