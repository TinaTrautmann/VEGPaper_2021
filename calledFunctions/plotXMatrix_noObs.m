function [sname] = plotXMatrix_noObs(dataNames, dataXX, lat, lon, colLabel, spth)
% plots matrix plot of correlationof XX experiments with observations + the scatter
% of the input data (calculates the gridwise correlation along ntix)
% for different variables
%
% returns the save name of the figure
%
% Input:
%     dataNames     = experiment names, in the order dataObs, data1, data2
%     dataUnc       = structure of observational uncertainties, fields = variables with timeseries for each grid; size(npix,ntix)
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
% colLims = {[0.4 1];[0.75 1]}; % for T/ET & Qslow/Q
% colLimDiff  = [-0.2 0.2];

% colLims = {[0 3];[0 2];[0 1]}; %for Q components
% colLimDiff  = [-1 1];

% colLims = {[0 3];[0 2];[0 1];[0 0.5]}; %for ET components
% colLimDiff  = [-0.5 0.5];

colLims     = colLabel.colLim;
colLimDiff  = colLabel.colLimDiff;

col         = othercolor('RdBu4',100);
colDiff     = othercolor('PuOr11',201);

pix_a = AreaGridLatLon(lat,lon, [1,1]);
pix_a = pix_a(:,1);

% show only the significant changes >0.05 (5%)
xx_sig = 0.05;
colDiff_sig = colDiff;
colDiff_sig(size(colDiff,1)/2-1:size(colDiff,1)/2+1,:) = [1,1,1;1,1,1;1,1,1];

expNames = dataNames;
dataNames = strrep(dataNames,'_','-');
dataNames2 = strrep(dataNames,'E-B-bL-RD4', 'VEG'); 



varNames = fieldnames(dataXX{1});
units = '-';

dN      = numel(dataXX);
nrows   = numel(dataXX);
ncols   = numel(dataXX);

% preps maps
land        = shaperead('landareas', 'UseGeoCoords', true);
rivers      = shaperead('worldrivers', 'UseGeoCoords', true);
geoRaRef    = georasterref('RasterSize', [180 360], 'RasterInterpretation', 'cells',  ...
    'LatitudeLimits', [-90 90], 'LongitudeLimits', [-180 180]);
Z1           = NaN(180, 360); % prepare mapgrid

%% Loop over variables
for vN = 1:numel(varNames)
    varN = varNames{vN};
    for dn = 1:dN
        expN = expNames{dn};
        tmpData = dataXX{dn};
        met1.(varN).(expN)    = nanmean(tmpData.(varN),2);
    end

end

%% Figure
    
    % loop over variables to plot
    for vN = 1:numel(varNames)
        varN  = varNames{vN};
        varN2 = strrep(varN,'_', ' / ');     
        
        colLim = colLims{vN};
        sp  = figure; set(gcf, 'Position', [5 5 ncols*5 nrows*4]);
        ha  = tight_subplot(nrows,ncols,[.02 .02],[.07 .05],[.05 .02]);
        % loop over rows (=variables), cols are added individually
        cnt = 1;
        for rr=1:nrows
            
            for cc=1:ncols
                
                if rr==cc % map of dataNames{rr}
                    
                    data        = met1.(varN).(expNames{cc});
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
                    
                    colormap(ha(cnt),col);
                    caxis(ha(cnt),colLim);
                    
                    t   =   title([ dataNames2{rr} ' (avg = ' num2str(round(nanmeanArea(data,pix_a),2)) ')'],'Fontsize', 7);
                    
                elseif rr<cc %difference map (rr-cc)
                    data        = met1.(varN).(expNames{rr})-met1.(varN).(expNames{cc});
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
                    
                    colormap(ha(cnt),colDiff);
                    caxis(ha(cnt),colLimDiff);
                    
                    title([dataNames2{rr} ' - ' dataNames2{cc}],'Fontsize', 7)
                    
                elseif rr>cc % histogram of difference cc-rr [OR: scatter of (cc,rr)]
                    axes(ha(cnt));
                    %histogram of metric difference
%                     c3 = met1.(varN).(expNames{cc})-met1.(varN).(expNames{rr});
%                     h=histogram(c3(:));
%                     set(ha(cnt), 'XLim', [prctile(c3, 5) prctile(c3, 95)], 'Fontsize', 6)
%                     h.FaceColor = rgb('Gray');
%                     hold on, plot([0 0], ha(cnt).YLim, '-', 'Color', rgb('DarkRed'), 'LineWidth', 1.5)
%                     ha(cnt).Position(1) = ha(cnt).Position(1) + 0.01;
%                     ha(cnt).Position(2) = ha(cnt).Position(2) + 0.03;
%                     ha(cnt).Position(3) = 0.15;
%                     ha(cnt).Position(4) = 0.13;

                    h=plot(met1.(varN).(expNames{cc})(:),met1.(varN).(expNames{rr})(:),'ok', 'MarkerSize', 3);
                    set(ha(cnt), 'XLim', [-0.02 1.02], 'YLim', [-0.02 1.02], 'Fontsize', 6)
                    hold on, plot([-0.02 1.02],[-0.02 1.02], '-', 'Color', rgb('DarkRed'), 'LineWidth', 1)
                    ha(cnt).Position(1) = ha(cnt).Position(1) + 0.04;
                    ha(cnt).Position(2) = ha(cnt).Position(2) + 0.08;
                    ha(cnt).Position(3) = 0.35;
                    ha(cnt).Position(4) = 0.35;
                    xlabel(dataNames2{cc}, 'Fontsize', 6),ylabel(dataNames2{rr}, 'Fontsize', 6)
                    
                end
                
                % add the colorbars
                if rr==1
                    if cc==1 %colorbar for correlation
                        cb  =   colorbar;
                        cb1 =   ylabel(cb, units);
                        set(cb,'FontSize',6, 'Position', [0.05    0.06    0.45    0.016],...
                            'Orientation', 'horizontal');
                        set(cb1,'Position',[(colLim(2)+colLim(1))/2 -2 0],'Fontweight', 'b');
                        
                    elseif cc==2 %colorbar for difference
                        cb  =   colorbar;
                        cb1 =   ylabel(cb, 'difference');
                        set(cb,'FontSize',6, 'Position', [0.525    0.06    0.45    0.016],...
                            'Orientation', 'horizontal');
                        set(cb1,'Position',[(colLimDiff(2)+colLimDiff(1))/2 -2 0],'Fontweight', 'b');
                    end
                end
                cnt = cnt+1;
            end
        end
        
        
        
        
        %     %add the overall title name
        an = annotation('textbox',[0 .5 1 .5],'String', [ strrep(varN2,'_','-') ] ,'LineStyle','none', 'Fontsize', 10, 'Fontweight', 'b','fitboxtotext', 'on', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center');
        
        % savename
        sname = ['Comparison Maps - ' varN ' - ' strjoin(dataNames2, ' vs ') ];
        
        if ~isempty(spth)
            print(gcf,[spth sname '_mean.png'],'-dpng','-r300');
        end
        
        close all
    end

end




