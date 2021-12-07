function [sname] = plotPieCR(expNames, varNames, CR_Avg,colVars, spth)
aggTimes    = {'msc','iav'};



ncols = numel(aggTimes);
nrows = numel(expNames);
cnt = 1;
sp  = figure; set(gcf, 'Position', [5 5 ncols*5 nrows*4]);
ha  = tight_subplot(nrows,ncols,[.02 .02],[.1 .1],[.1 .1]);

for eN=1:numel(expNames)
    expN = expNames{eN};
    for aN = 1:numel(aggTimes)
        agg = aggTimes{aN};
        
        dataExp = [];
        for vN=1:numel(varNames)
            varN = varNames{vN};
            tmp = CR_Avg.(expN).(agg).(varN);
            if tmp<=0; tmp=10^-3;end
            dataExp = [dataExp, tmp];
        end
        
        axes(ha(cnt))
        p=pie(dataExp);
        set(findobj(p,'type','text'),'fontsize',7)
        colormap(colVars);
        
        if aN==1
            tmpN = expN;
            tmpN       = strrep(tmpN,'_','-');
            if strcmp(tmpN,'E-B-bL-RD4')
                tmpN = 'VEG';
            end
            t1 = text(-0.2,0.55, strrep(tmpN,'_','-'), 'Units','Normalized', 'Fontsize', 10, 'Fontweight', 'b', 'HorizontalAlignment', 'center');
            t1.Rotation = 90;
        end
        
        if eN==1
            if strcmp(agg,'msc')
                tmpT = 'Mean Seasonal Cycle';
            elseif strcmp(agg,'iav')
                tmpT = 'Inter-Annual Variability';
            end

            t2 = text(1.1,1.2, tmpT, 'Units','Normalized', 'Fontsize', 10, 'Fontweight', 'b', 'HorizontalAlignment', 'right');
        end
        
        cnt = cnt+1;
    end
    
    
end

varNamesTmp = strrep(varNames,'I','I-')
varNamesTmp = strrep(varNamesTmp,'I','I-');
varNamesTmp = strrep(varNamesTmp,'wGW','wDeep');
varNamesTmp = strrep(varNamesTmp,'wSurf','wSlow');
l = legend(varNamesTmp, 'Position', [.35 .01 .3 .05], 'Orientation', 'Horizontal', 'Box', 'off');

% savename
sname = ['SpatialAvg_CR - ' strjoin(varNames, ' vs ') ];

if ~isempty(spth)
    print(gcf,[spth sname '.png'],'-dpng','-r300');
end

end