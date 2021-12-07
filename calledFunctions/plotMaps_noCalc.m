function [sname] = plotMaps_noCalc(dataNames,varNames, dataIn, lat, lon, colLabel, spth)
% plots maps of 2 experiments  + tboxplot of distributions difference
% of the input data (calculates the gridwise correlation along ntix)
% for different variables
%
% returns the save name of the figure
%
% Input:
%     dataNames     = experiment names, in the order data1, data2
%     varNames      = variable names/ fieldnames to plot for each (=2)
%     experiment
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


% preps
if numel(varNames)==5
    col_box     = [rgb('Blue');rgb('DarkRed');rgb('SaddleBrown');rgb('Green');rgb('Purple')]
elseif numel(varNames)==4
    col_box     = [rgb('Blue');rgb('SaddleBrown');rgb('Green');rgb('Purple')]
elseif numel(varNames)==3
    col_box     = [rgb('Blue');rgb('SaddleBrown');rgb('Purple')]
end
col         = othercolor('Blues9',100);
colDiff     = flipud(othercolor('PuOr11',201));
units = [];
pix_a = AreaGridLatLon(lat,lon,[1,1]);
pix_a = pix_a(:,1);
dataNames_org   = dataNames;
dataNames       = strrep(dataNames,'_','-');
dataNames = strrep(dataNames,'E-B-bL-RD4', 'VEG');

% idx             = strcmp(dataNames,'E-B-bL-RD4');
% dataNames{idx==1} = 'VEG';

varNames_org = varNames;
varNames = strrep(varNames,'I','I-');
varNames = strrep(varNames,'wGW','wDeep');
varNames = strrep(varNames,'wSurf','wSlow');

nrows   = numel(dataNames);
ncols   = numel(varNames)+1;

% preps maps
land        = shaperead('landareas', 'UseGeoCoords', true);
rivers      = shaperead('worldrivers', 'UseGeoCoords', true);
geoRaRef    = georasterref('RasterSize', [180 360], 'RasterInterpretation', 'cells',  ...
    'LatitudeLimits', [-90 90], 'LongitudeLimits', [-180 180]);
Z1           = NaN(180, 360); % prepare mapgrid




%% Figure
% loop over msc/ to plot
metNames = {'msc','iav'};

for cn = 1:numel(metNames)
    mN       = metNames{cn};
    dataCell = cell(1,numel(varNames));
    
    sp  = figure; set(gcf, 'Position', [5 5 ncols*5 nrows*4]);
    ha  = tight_subplot(nrows,ncols,[.02 .02],[.05 .05],[.05 .02]);
    % loop over rows (=experiments), cols are added individually
    cnt = 1;
    for rr=1:nrows
        expN_org = dataNames_org{rr};
        expN = dataNames{rr};
        
        for cc=1:ncols
            axes(ha(cnt))
            
            if cc==numel(varNames)+1 %(boxplot of the metrics)
                
                %  %boxplot metric values
                aboxplot(dataCell ,'labels', varNames, 'colormap', col_box, 'widths', 0.4)
                set(ha(cnt), 'XTickLabel', [], 'FontSize', 6, 'YTickLabelMod', 'auto')
                hold on, plot(ha(cnt).XLim,[0 0],'k-');
%                 ha(cnt).Position = [0.79  1-rr*0.3 0.18 0.2]
%                 ha(cnt).Position(1) = 0.82;
%                 ha(cnt).Position(2) = ha(cnt).Position(2) + 0.09
%                 ha(cnt).Position(3) = 0.16;
%                 ha(cnt).Position(4) = 0.30;
%                 set(ha(cnt), 'YLim', [-1 1])
                if rr==1
                    title(['Distribution'], 'Fontsize', 7);
                end

                for iN =1:numel(dataCell)
                    tmpI = nanmeanArea(dataCell{iN},pix_a);
                    t1 = text(0.82,0.55+0.1*iN, ['$\bar{s}$' '=' num2str(round(tmpI,2)) ],'Units','Normalized', 'Fontsize', 6, 'Color', col_box(iN,:), 'HorizontalAlignment', 'left', 'Interpreter', 'Latex');
                end
                
                if rr==nrows
                    %                     l = legend(dataNames(2:3), 'Box', 'off', 'AutoUpdate', 'off');
                    l = legend(varNames, 'Box', 'off', 'AutoUpdate', 'off', 'Orientation', 'Horizontal');
                    l.Position = [0.74 0.03 0.19 0.01];
                end
                
                
            else
                
                dataCell{cc} = dataIn.(expN_org).(mN).(varNames_org{cc});
                data         =  dataIn.(expN_org).(mN).(varNames_org{cc});
                
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
                
                
                t   =   title([ num2str(cc) ') ' varNames{cc}],'Fontsize', 7);
                colLim = [0 1];
                colormap(ha(cnt),col);
                caxis(ha(cnt),colLim);
                if rr==nrows
                    if cc==1
                        cb  =   colorbar;
                        cb1 =   ylabel(cb, units);
                        set(cb,'FontSize',6, 'Position', [0.05 0.05 0.55 0.015],...
                            'Orientation', 'horizontal');
                        set(cb1,'Position',[(colLim(2)+colLim(1))/2 -2 0]);
                    end
                end
            end
            
            %add the experiment name in the first column
            if cc==1
                t1 = text(-0.1,0.55, [expN], 'Units','Normalized', 'Fontsize', 10, 'Fontweight', 'b', 'HorizontalAlignment', 'center');
                t1.Rotation = 90;
            end
            cnt = cnt+1;
        end
        
        
        
    end
    
    %     %add the overall title name
    if strcmp(mN,'msc')
    tmpT = 'Mean Seasonal Cycle';
    elseif strcmp(mN,'iav')
    tmpT = 'Inter-Annual Variability';
    end
    an = annotation('textbox',[0 .5 1 .5],'String', ['TWS Composition - ' tmpT ] ,'LineStyle','none', 'Fontsize', 10, 'Fontweight', 'b','fitboxtotext', 'on', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center');
    
    %' - ' strjoin(varNames, ' vs ')
    
    
    % savename
    sname = ['Comparison ' mN ' Maps - ' strjoin(varNames, ' vs ') ];
    
    if ~isempty(spth)
        print(gcf,[spth sname '.png'],'-dpng','-r300');
    end
    
end

end


