function [sname] = plotDiffMaps_InputPercDiff(dataNames, varNames, dataIn, lat, lon, colLabel, spth)
% plots % difference maps of 2 experiments  + histogram of the difference
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
col         = othercolor('Blues9',100);
colDiff     = flipud(othercolor('BrBG10',401));
colDiff     = colDiff(200:end,:);


units = [];

dataNames_org = dataNames;
dataNames = strrep(dataNames,'_','-');
dataNames = strrep(dataNames, 'E-B-bL-RD4', 'VEG');

nrows   = numel(varNames);
ncols   = 4;

% preps maps
land        = shaperead('landareas', 'UseGeoCoords', true);
rivers      = shaperead('worldrivers', 'UseGeoCoords', true);
geoRaRef    = georasterref('RasterSize', [180 360], 'RasterInterpretation', 'cells',  ...
    'LatitudeLimits', [-90 90], 'LongitudeLimits', [-180 180]);
Z1           = NaN(180, 360); % prepare mapgrid




%% Figure
    sp  = figure; set(gcf, 'Position', [5 5 ncols*5 nrows*4]);
    ha  = tight_subplot(nrows,ncols,[.02 .02],[.05 .05],[.05 .02]);
    % loop over rows (=variables), cols are added individually
    cnt = 1;
    for rr=1:nrows
        
        varN = varNames{rr};
        
        c1 = dataIn.(dataNames_org{1}).(varN);
        c2 = dataIn.(dataNames_org{1}).(varN) - dataIn.(dataNames_org{2}).(varN);
        c3 = (c2./c1)*100;
        

        for cc=1:ncols
            
            if cc==4 %(boxplot of the metrics)
                axes(ha(cnt))

               %  histogram of difference
%                 h=histogram(c3)
                d2 = histcounts(c3);
                b = bar(d2, 'facecolor', 'flat');
                colHist     = flipud(othercolor('BrBG10',2*length(b.XData)));
                colHist     = colHist(length(b.XData)+1:end,:);

                b.CData = colHist;
                b.EdgeColor = 'none';
                set(ha(cnt), 'FontSize', 6, 'YTickLabelMod', 'auto', 'box', 'on')
                ha(cnt).Position(1) = 0.79;
                ha(cnt).Position(2) = 0.13;
                ha(cnt).Position(3) = 0.18;
                ha(cnt).Position(4) = 1/nrows-0.36;
                title(['% Difference'], 'Fontsize', 7);   
                cl=xlabel('%')
                cl.Position(2) = -80;
                ylabel('# of grid cells')
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
                
                if cc==3 % = percent difference data1-data2
                    t   =   title('% Difference','Fontsize', 7);
                    caxis([0 30])
                    colormap(ha(cnt),colDiff);
                    cb  =   colorbar;
                    set(cb,'FontSize',6,...
                        'Orientation', 'horizontal', 'Location', 'South', 'AxisLocation', 'out');
                    cb.Position(2) = cb.Position(2)-0.13;
                    cb.Position(3) = 0.2;
                    cl = ylabel(cb, '%');
                    cl.Position(2) = -1;

                else % cols 1 & 2
                    if cc==1
                        t   =   title([ dataNames{cc}],'Fontsize', 7);
                        colormap(ha(cnt),col);

                    elseif cc==2
                        t   =   title([ 'Abs Difference ' dataNames{cc-1} ' - ' dataNames{cc}],'Fontsize', 7);
                        colormap(ha(cnt),colDiff);
                    end
  
                    cb  =   colorbar;
                    set(cb,'FontSize',6,...
                        'Orientation', 'horizontal', 'Location', 'South', 'AxisLocation', 'out');
                    cb.Position(2) = cb.Position(2)-0.13;
                    cb.Position(3) = 0.2;
                    cl = ylabel(cb, 'mm');
                    cl.Position(2) = -1;
                end
                
                %add the experiment name in the first column
                if cc==1
                    varN2 = 'mean daily ET'
                    t1 = text(-0.1,0.55, [strrep(varN2,'_','-')], 'Units','Normalized', 'Fontsize', 10, 'Fontweight', 'b', 'HorizontalAlignment', 'center');
                    t1.Rotation = 90;
                end
                
            end
            
            cnt = cnt+1;
            
        end
        
        %     %add the overall title name
%         an = annotation('textbox',[0 .5 1 .5],'String', ['Comparison of TWS Contributions ' strrep(mN,'_','-') ' - ' varNames{1} ' vs ' varNames{2}] ,'LineStyle','none', 'Fontsize', 10, 'Fontweight', 'b','fitboxtotext', 'on', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center');
        
    end
    
    
    % savename
    sname = ['Comparison Maps PercDiff - ' strjoin(dataNames_org, ' vs ') ' - ' strjoin(varNames, '_')];
    
    if ~isempty(spth)
        print(gcf,[spth sname '.png'],'-dpng','-r300');
    end
    

end


