%% Plot MSC of obs & 2 simulations for 3 variables globally & for different zones
function [sname, T_corr, sp] = plotMSCvarsZones_xExp_noObs(dataNames, dataXX, zoneNames, zoneIdx, pix_a, unitNames, col)
% plots the global and zonal MSC of x experiments against observations
% for different variables
%
% returns the save name of the figure, the figure handle and a
% table-structure with the correlation of the 1-2 experiments with
% observations for all zones and variables
%
% Input:
%     dataNames     = experiment names, in the order dataObs, data1, data2
%     dataXX        = cellarray with structures of each experiment, fields = variables with MSC for each grid; size(npix,12)
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
    col         = [rgb('Blue'); rgb('Green')];
end


varNames = fieldnames(dataXX{1});


dataNames  = strrep(dataNames,'_','-');
dataNames2 = strrep(dataNames,'E-B-bL-RD4', 'VEG'); 
zoneNames = strrep(zoneNames,'_','-');

eN = numel(dataXX);


% varNames  = {'ET', 'Q'}
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
end

sp  = figure; 
if nrows==1
    set(gcf, 'Position', [5 1 ncols*5 nrows*5+2]);
    ha  = tight_subplot(nrows,ncols,[.05 .02],[.15 .15],[.05 .02]);
else
    set(gcf, 'Position', [5 1 ncols*5 nrows*5]);
    ha  = tight_subplot(nrows,ncols,[.05 .02],[.08 .08],[.05 .02]);
end
% loop over rows (=variables), then cols (=global + zones)
cnt = 1;
for rr=1:nrows
    varN  = varNames{rr};
    varN2 = strrep(varN,'_', ' / ');
    unitN = unitNames{rr};
    for dN = 1:eN
        tmp = dataXX{dN};
        eval(char(['d' num2str(dN) '= tmp.(varN);']));
    end

    for cc=1:ncols
        if cc==1
            idxZ = 1:1:size(pix_a,1);
        else
            idxZ = find(zoneIdx==cc-1);
        end
        axes(ha(cnt));
        
        for dN = 1:eN
            eval(char(['dXX = d' num2str(dN) ';']))
            plot(xSeason,nanmeanArea(dXX(idxZ,:),pix_a(idxZ)), 'LineStyle', '-', 'Color', col(dN,:)), hold on
        end
        plot([0.5 12.5],[0 0], '-', 'color', [.5 .5 .5])
        ylim([0 1])
        if rr==1; title(zNames{cc}); end
        if cc==1; yl= ylabel([varN2 ' [' unitN ']'], 'Fontweight', 'b'); yl.Position(1) = -1.5; end
        ha(cnt).XAxis.MinorTickValues =  [1:1:12];
        cnt = cnt+1;
          
    end
end
% set axis
set(ha(:), 'Xlim', [0.5 12.5], 'XTick', [3:3:12], 'GridAlpha', 0.25, 'XMinorGrid', 'on', 'MinorGridAlpha', 0.05, 'MinorGridLineStyle', '-',...
    'YGrid', 'on', 'Fontsize', 7)
set(ha(1:end-ncols), 'XTickLabel', [])

for tmp=1:length(ha),ha(tmp).YLabel.FontSize = 10;end

% legend & overall title
tmpT = dataNames2;

if nrows==1
    l  = legend(tmpT, 'Position', [.44 0 .12 .07], 'box', 'off', 'Orientation', 'Horizontal','Fontsize', 7);
else                    
    l  = legend(tmpT, 'Position', [.44 -.02 .12 .07], 'box', 'off', 'Orientation', 'Horizontal','Fontsize', 7);
end
an = annotation('textbox',[0 .5 1 .5],'String', ['Comparison of the Mean Seasonal Cycle - ' strjoin(tmpT, ' vs ') ] ,'LineStyle','none', 'Fontsize', 10, 'Fontweight', 'b','fitboxtotext', 'on', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center');

% savename
sname = ['Comparison of the Mean Seasonal Cycle - ' strjoin(tmpT, ' vs ') ];

end