%% Plot Parameter Uncertainty
%% set save pth
spth = [pwd '/plots/ParameterUncertainty/'];
if ~exist(spth, 'dir'), mkdir(spth),end


%% define models
expNames    = {'B', 'E_B_bL_RD4', 'PFT'};
rDates      = {'20200629', '20200629', '20210917'};

space_run   =  'validCali';

%% load the stuff
for eN = 1:numel(expNames)
    expName = expNames{eN};
    rDate   = rDates{eN};
    
    load(['data/model_output/' space_run '_' expName '_' space_run '_' rDate '/optimization/Results_ParamUncertainty_' expName '.mat'])
    
    %plot the correlation matrix
    pa_name     = Tpara_se.Properties.RowNames;
    nP          = numel(pa_name);
    paraName    = cell(1,nP);
    for pn=1:nP
        tmp             = strsplit(pa_name{pn},'.');
        paraName{pn}    = strrep(tmp{end},'_','-');
    end
    
    cols = othercolor('RdBu11');    
    figure, imagesc(cormatDiffSig);
    set(gca, 'YTick', 1:nP, 'YTickLabel', paraName, 'XTick', 1:nP, 'XTickLabel', paraName, 'XAxisLocation', 'top', 'XTickLabelRotation', 90, 'fontsize', 8)
    set(gcf, 'Position', [5 5 15 13])
    set(gca, 'Position', [.16 .02 .7 .8])
    cb      = colorbar;
    colormap(cols)
    set(cb, 'Limits', [-1 1])
    ylabel(cb, 'correlation')
    t= title(['Param Correlation - ' strrep(expName,'_','-') ], 'fontsize', 10, 'fontweight', 'b', 'Position', [5.7 -1.5 .8]);
    print(gcf,[spth ['Param Correlation DiffSig - ' char(expName) ] '.png'],'-dpng','-r300');
    
    
    % -- only show correlations > 0.5
    cols = othercolor('RdBu11',20);
    cols(6:15,:) = 1;
    xx=0.5;
    tmpC = cols;
    tmpC(5,:) = 1;
    
    tmpDiffSig = zeros(size(cormatDiffSig));
    tmpDiffSig(abs(cormatDiffSig)>=0.5) = cormatDiffSig(abs(cormatDiffSig)>=0.5);
    figure, imagesc(tmpDiffSig)
    set(gca, 'YTick', 1:nP, 'YTickLabel', paraName, 'XTick', 1:nP, 'XTickLabel', paraName, 'XAxisLocation', 'top', 'XTickLabelRotation', 90, 'fontsize', 8)
    set(gcf, 'Position', [5 5 15 13])
    set(gca, 'Position', [.16 .02 .7 .8])
    cb      = colorbar;
    % tmpC(4,:) = 1;
    colormap(tmpC)
    set(cb, 'Limits', [-1 1])
    ylabel(cb, 'correlation')
    t= title(['Param Correlation >= 0.5 - ' strrep(expName,'_','-') ], 'fontsize', 10, 'fontweight', 'b', 'Position', [5.7 -1.5 .8]);
    print(gcf,[spth ['Param Correlation DiffSig_05 - ' char(expName) ] '.png'],'-dpng','-r300');
    
    % show the table
    Tpara_se
    
end
