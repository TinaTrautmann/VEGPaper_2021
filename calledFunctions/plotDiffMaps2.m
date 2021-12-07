function [sname] = plotDiffMaps2(dataNames, dataObs, data1, data2, lat, lon, colLabel, spth)
% plots correlation maps of 1-2 experiments with observations + the scatter
% of the input data (calculates the gridwise correlation along ntix)
% for different variables
%
% returns the save name of the figure
%
% Input:
%     dataNames     = experiment names, in the order dataObs, data1, data2
%     dataObs       = structure of observations, fields = variables with timeseries for each grid; size(npix,ntix)
%     data1         = structure of experiment1, fields = variables with timeseries for each grid; size(npix,ntix)
%     data2         = structure of experiment2, fields = variables with timeseries for each grid; size(npix,ntix)
%     lat           = latitude per grid;
%     lon           = longitude per grid;
%     colLabel.size = size of the points in the grid if only a global subset
%     should be plotted
%     colLimScatter = x and y limits for the scatter plot; per variable,
%     e.g. colLimScatter = {[0 5], [0 5], [0 5]}; %per variable

%% plotDiffMaps
set(0,'DefaultLegendAutoUpdate','off')

% plot global maps of 2 experiment, their difference and a scatterplot


% %% EXAMPLE Input = original (monthly) data, 3 variables
% load('data/input/globalBaseline_validCali/globalBaseline_Constraints_1deg.mat', 'lat', 'lon')
% dataNames = {'obs', 'exp1', 'exp2'};
% dataObs.A = randi(5,length(lat),12*5);
% data1.A   = dataObs.A ./2 + randi(2,length(lat),12*5);
% data2.A   = max(0,dataObs.A - randi(2,length(lat),12*5));
%
% dataObs.B = randi(5,length(lat),12*5);
% data1.B   = dataObs.B ./2 + randi(2,length(lat),12*5);
% data2.B   = max(0,dataObs.B - randi(2,length(lat),12*5));
%
% dataObs.C = randi(5,length(lat),12*5);
% data1.C   = dataObs.C ./2 + randi(2,length(lat),12*5);
% data2.C   = max(0,dataObs.C - randi(2,length(lat),12*5));
%
% colLabel.size = [];
% colLimScatter = {[0 5], [0 5], [0 5]};


% preps
col         = othercolor('RdBu4',100);
colDiff     = othercolor('PuOr11',201);
% show only the significant changes >0.05 (5%)
xx_sig = 0.05;
colDiff_sig = colDiff;
colDiff_sig(size(colDiff,1)/2-1:size(colDiff,1)/2+1,:) = [1,1,1;1,1,1;1,1,1];


dataNames = strrep(dataNames,'_','-');
dataNames = strrep(dataNames,'E-B-bL-RD4', 'VEG');
varNames    = fieldnames(dataObs);

nrows   = numel(varNames);
ncols   = 4;

% preps maps
land        = shaperead('landareas', 'UseGeoCoords', true);
rivers      = shaperead('worldrivers', 'UseGeoCoords', true);
geoRaRef    = georasterref('RasterSize', [180 360], 'RasterInterpretation', 'cells',  ...
    'LatitudeLimits', [-90 90], 'LongitudeLimits', [-180 180]);
Z1           = NaN(180, 360); % prepare mapgrid


%% calculate gridwise correlation (can add other metrics in columns and make them a table
for rr=1:numel(varNames)
    varN  = varNames{rr};
    % correlation
        met1.corr.(varN) = NaN(length(lat),1);
        met2.corr.(varN) = NaN(length(lat),1);
    
    % KGE
    met1.KGE.(varN) = NaN(length(lat),1);
    met1.corr_P.(varN) = NaN(length(lat),1);
    met1.var_err.(varN) = NaN(length(lat),1);
    met1.bias_err.(varN) = NaN(length(lat),1);
    
    met2.KGE.(varN) = NaN(length(lat),1);
    met2.corr_P.(varN) = NaN(length(lat),1);
    met2.var_err.(varN) = NaN(length(lat),1);
    met2.bias_err.(varN) = NaN(length(lat),1);
%     % MEF
    
    % loop over grids
    for np=1:length(lat)
        
        %Correlation
        met1.corr.(varN)(np,1) = corr(dataObs.(varN)(np,:)',data1.(varN)(np,:)','rows', 'complete');
        met2.corr.(varN)(np,1) = corr(dataObs.(varN)(np,:)',data2.(varN)(np,:)','rows', 'complete');
        
        % KGE
        [KGE, r_p, alpha, beta] = calcKGE(dataObs.(varN)(np,:)',data1.(varN)(np,:)');
        
        met1.KGE.(varN)(np,1) = KGE;
        met1.corr_P.(varN)(np,1) = r_p;
        met1.var_err.(varN)(np,1) = (alpha-1)^2;
        met1.bias_err.(varN)(np,1) = (beta-1)^2;
        
        
        [KGE, r_p, alpha, beta] = calcKGE(dataObs.(varN)(np,:)', data2.(varN)(np,:)');
        met2.KGE.(varN)(np,1) = KGE;
        met2.corr_P.(varN)(np,1) = r_p;
        met2.var_err.(varN)(np,1) = (alpha-1)^2;
        met2.bias_err.(varN)(np,1) = (beta-1)^2;
        
        % MEF
        v_obs   = find(~isnan(dataObs.(varN)(np,:)));

        NSE = calcMEF(dataObs.(varN)(np,v_obs)',data1.(varN)(np,v_obs)', ones(size(data1.(varN)(np,v_obs)')));
        met1.MEF.(varN)(np,1) = 1-NSE;
        
        NSE = calcMEF(dataObs.(varN)(np,v_obs)',data2.(varN)(np,v_obs)', ones(size(data2.(varN)(np,v_obs)')));
        met2.MEF.(varN)(np,1) = 1-NSE;
    end
end



%% Figure
% loop over metrics to plot
metNames = fieldnames(met1);
metNames = {'corr_P'};
for cn = 1:numel(metNames)
    mN      = metNames{cn};
    
    
    sp  = figure; set(gcf, 'Position', [5 5 ncols*5 nrows*4]);
    ha  = tight_subplot(nrows,ncols,[.02 .02],[.05 .05],[.05 .02]);
    % loop over rows (=variables), cols are added individually
    cnt = 1;
    for rr=1:nrows
        varN  = varNames{rr};
        
        c1   = met1.(mN).(varN);
        c2   = met2.(mN).(varN);
        c3   = c1-c2;
        units = 'Pearson Correlation';
%         units   = strrep(mN,'_','-');
        
        % define color limits
        if strcmp(mN,'KGE') || strcmp(mN,'corr_P') || strcmp(mN,'MEF') || strcmp(mN,'corr') ==1
            colLim      = [0 1];
            colLimDiff  = [-0.2 0.2];
        else
            colLim      = [0 max(prctile(c1,90),prctile(c2,90))];
            colLimDiff  = [-0.5 0.5];
        end
        
        for cc=1:ncols
            
            if cc==4 %(boxplot of the metrics)
                axes(ha(cnt))
                %                 %boxplot metric values
                %                 aboxplot({c1,c2} ,'labels', dataNames(2:3), 'colormap', [rgb('DarkBlue');rgb('Green')], 'widths', 0.4)
                %                 set(ha(cnt),'Position', [0.79  1-rr*0.3 0.18 0.2], 'XTickLabel', [], 'FontSize', 6, 'YTickLabelMod', 'auto')
                %                 if strcmp(mN,'KGE') || strcmp(mN,'corr_P')==1
                %                     set(ha(cnt),'YTick',[0:0.2:1], 'YTickLabel', [0:0.2:1], 'YLim', [0 1.1])
                %                 end
                
%                 %boxplot metric difference
%                 aboxplot({c3} ,'labels', {'difference'}, 'colormap', [rgb('Gray')], 'widths', 0.4)
%                 set(ha(cnt),'Position', [0.79  1-rr*0.3 0.18 0.2], 'XTickLabel', [], 'FontSize', 6, 'YTickLabelMod', 'auto')
%                 set(ha(cnt), 'YLim', [prctile(c3, 5) prctile(c3, 95)])
                
                %histogram of metric difference
                h=histogram(c3);
                set(ha(cnt),'Position', [0.79  1-rr*0.3 0.18 0.2], 'FontSize', 6, 'YTickLabelMod', 'auto')
                set(ha(cnt), 'XLim', [prctile(c3, 5) prctile(c3, 95)])
                h.FaceColor = rgb('Gray');
                hold on, plot([0 0], ha(cnt).YLim, '-', 'Color', rgb('DarkRed'), 'LineWidth', 1.5)
                %             if rr==1
                ylabel('# of grid cells')
                title(['Distribution of Difference' ], 'Fontsize', 7);
                %             end
%                 if rr==nrows
%                     %                     l = legend(dataNames(2:3), 'Box', 'off', 'AutoUpdate', 'off');
%                     l = legend({'1) - 2)'}, 'Box', 'off', 'AutoUpdate', 'off');
%                     l.Position = [0.79 0.058 0.19 0.02];
%                 end
            else
                data = eval(char(['c' num2str(cc)]));
                Z           = NaN(180, 360); % prepare mapgrid
                Z           = imbedm(lat, lon, data, Z, geoRaRef); % insert data in mapgrid
                Z           = [Z ones(180,1)];  % add column for correct plotting of last column
                Z           = [Z; ones(1,361)]; % add row for correct plotting of last row
                
                axes(ha(cnt));
                worldmap('World');
                ha(cnt).Layer = 'top';
                setm(gca,'ParallelLabel', 'off','MeridianLabel','off', 'FLineWidth', 0.5);
                geoshow(ha(cnt), land, 'FaceColor', [0.85 0.85 0.85])
                hold on
                if isfield(colLabel,'size')
                    geoshow(ha(cnt), rivers, 'Color', [.25 .25 .25]);
                    sm = scatterm(lat, lon, colLabel.size, data, 'filled');
                else
                    sm = surfm([-90 90], [-180 180], Z);
                    geoshow(ha(cnt), rivers, 'Color', [.25 .25 .25]);
                end
                
                if cc==3 %= difference data1-data2
                    
                    %                     data_sig = data;
                    %                     data_sig(data_sig>-xx_sig&data_sig<xx_sig) = 0;
                    %
                    %                     if isfield(colLabel,'size')
                    %                         geoshow(ha(cnt), rivers, 'Color', [.25 .25 .25]);
                    %                         sm = scatterm(lat, lon, colLabel.size, data_sig, 'filled');
                    %                     else
                    %                         sm = surfm([-90 90], [-180 180], Z);
                    %                         geoshow(ha(cnt), rivers, 'Color', [.25 .25 .25]);
                    %                     end
                    %                     colormap(ha(cnt),colDiff_sig);
                    
                    colormap(ha(cnt),colDiff);
                    
                    
                    %                 if rr==1
                    t   =   title('1) - 2)','Fontsize', 7);
                    %                 end
                    
                    
                    if strcmp(mN,'KGE') || strcmp(mN,'corr_P') || strcmp(mN,'MEF') || strcmp(mN,'corr')==1
                        caxis(ha(cnt),colLimDiff);
                        if rr==nrows
                            cb  =   colorbar;
                            cb1 =   ylabel(cb, 'difference');
                            set(cb,'FontSize',6, 'Position', [0.533 0.07 0.2 0.0175],...
                                'Orientation', 'horizontal');
                            set(cb1,'Position',[(colLimDiff(2)+colLimDiff(1))/2 -2 0]);
                        end
                    else
                        cb  =   colorbar;
                        caxis(ha(cnt),colLimDiff);
                        set(cb,'FontSize',6, 'Position', [0.533 1-rr*0.31 0.2 0.015],...
                            'Orientation', 'horizontal');
                        if rr==nrows
                            cb1 =   ylabel(cb, 'difference');
                        end
                    end
                    
                else
                    
                    tmpT = dataNames{cc+1};
                    if strcmp(tmpT,'E-B-bL-RD4')
                        tmpT = 'VEG';
                    end
                    %                 if rr==1
                    t   =   title([ num2str(cc) ') ' tmpT],'Fontsize', 7);
                    %                 end
                    if strcmp(mN,'KGE') || strcmp(mN,'corr_P') || strcmp(mN,'MEF') || strcmp(mN,'corr') ==1
                        colormap(ha(cnt),col);
                    else
                        colormap(ha(cnt),flipud(col));
                    end
                    caxis(ha(cnt),colLim);
                    if  strcmp(mN,'KGE') || strcmp(mN,'corr_P') || strcmp(mN,'MEF') || strcmp(mN,'corr') ==1
                        if rr==nrows
                            if cc==1
                                cb  =   colorbar;
                                cb1 =   ylabel(cb, units);
                                set(cb,'FontSize',6, 'Position', [0.05 0.07 0.45 0.015],...
                                    'Orientation', 'horizontal');
                                set(cb1,'Position',[(colLim(2)+colLim(1))/2 -2 0]);
                            end
                        end
                    else
                        cb  =   colorbar;
                        set(cb,'FontSize',6, 'Position', [0.06+0.236*(cc-1) 1-rr*0.31 0.2 0.015],...
                            'Orientation', 'horizontal');
                        if rr==nrows
                            cb1 =   ylabel(cb, units);
                        end
                    end
                end
                
                %add the variable name in the first column
                if cc==1
                    t1 = text(-0.1,0.55, [varN], 'Units','Normalized', 'Fontsize', 10, 'Fontweight', 'b', 'HorizontalAlignment', 'right');
                    t1.Rotation = 90;
                end
                
            end
            
            cnt = cnt+1;
            
        end
        

        %     %add the overall title name
%         an = annotation('textbox',[0 .5 1 .5],'String', ['Pearson Correlation - B vs VEG' ] ,'LineStyle','none', 'Fontsize', 10, 'Fontweight', 'b','fitboxtotext', 'on', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center');
        %strjoin(dataNames(2:3), ' vs ') 
    end
    
        an = annotation('textbox',[0 .5 1 .5],'String', ['Comparison of ' strrep(mN,'_','-')  ' - ' strjoin(dataNames(2:3), ' vs ') ] ,'LineStyle','none', 'Fontsize', 10, 'Fontweight', 'b','fitboxtotext', 'on', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center');
    
    % savename
    sname = ['Comparison ' mN ' Maps - ' strjoin(dataNames, ' vs ') ];
    
    if ~isempty(spth)
        print(gcf,[spth sname '.png'],'-dpng','-r300');
    end
    
end

end


