function [ sp ] = PlotMapGlobal_noInfo(name, addinfo, units, lat, lon, data, colLim, col, colLabel, cbar, histPlot)
% plots data onto a global Map also land area and large rivers
% returns the save name
%       name                - title of the Map
%       addinfo             - additional information
%       units               - units of colorbar
%       lat                 - latitude (pix,1)
%       lon                 - longitude (pix,1)
%       data                - data to plot (pix,1)
%       colLim              - limits of the colorbar, default = [min prctile(data,95)]
%       col                 - colors for the colorbar, default = jet; example: col=[.7,0,0; .5,.5,.5; 0,.7,0]; OR col=flipud(jet)
%       colLabel            - structure containing the yticks and y
%       colLabel.cticks     = [-0.75,0,0.75]
%       colLabel.clabels    = {'orig', 'no difference', 'ETconstAlpha'}
%       colLabel.size       = markersize of points -> uses scatterm instead of surfm to plot the data
%       cbar                = 1 to plot the colorbar, 0 to not plot it
%       histPlot            = 1 to plot the histogram of the data

%
land    = shaperead('landareas', 'UseGeoCoords', true);
rivers  = shaperead('worldrivers', 'UseGeoCoords', true);


geoRaRef    = georasterref('RasterSize', [180 360], 'RasterInterpretation', 'cells',  ...
    'LatitudeLimits', [-90 90], 'LongitudeLimits', [-180 180]);
Z           = NaN(180, 360); % prepare mapgrid
Z           = imbedm(lat, lon, data, Z, geoRaRef); % insert data in mapgrid
Z           = [Z ones(180,1)];  % add column for correct plotting of last column
Z           = [Z; ones(1,361)]; % add row for correct plotting of last row

if isempty(colLim)==1
    colLim  = [floor(min(data(:))) ceil(prctile(data(:),95))];
    if colLim(1) == colLim(2), colLim(2)=colLim(2)+2;end
end

if isempty(col)==1
    col     = jet;
end

if isempty(cbar)==1
    cbar     = 1;
end

if histPlot == 1
    sb1=subplot(1,2,1);
end
ax=worldmap('World'); 
% latlim  = [30 90];
% lonlim  = [-180 180];
% ax = worldmap([latlim lonlim]);

set(ax, 'layer', 'top','xcolor', 'w', 'ycolor', 'w'); %'Position', [0.1 0.12 0.8 0.75]
axis off; grid on; axis tight;
setm(gca,'ParallelLabel', 'off','MeridianLabel','off', 'FLineWidth', 0.5);
geoshow(ax, land, 'FaceColor', [0.85 0.85 0.85]);
hold on

if isfield(colLabel,'size')
    geoshow(ax, rivers, 'Color', [.25 .25 .25]);
    sm = scatterm(lat, lon, colLabel.size, data, 'filled');
else
    sm = surfm([-90 90], [-180 180], Z);
    geoshow(ax, rivers, 'Color', [.25 .25 .25]);
end

caxis(colLim);
colormap(col);
t   =   title(name);

if cbar == 1
    cb  =   colorbar;
    cb1 =   ylabel(cb, units);
    % edit colorbar
    set(cb,'FontSize',8, 'Position', [0.1 0.09 0.8 0.02],...
        'Orientation', 'horizontal');
    try  % set limits and ticknames if they are specified
        cticks  = colLabel.cticks;
        clabels = strrep(colLabel.clabels,'_','-');
        set(cb, 'YTick', cticks, 'YTickLabel', clabels);
    end
    set(cb1,'Position',[(colLim(2)+colLim(1))/2 -2 0]);
end

% title and heading
set(t,'Fontsize', 9, 'Fontweight', 'b', 'Position', [9 9700000 0 ]); %worldmap
annotation('textbox', 'String', addinfo, 'fontsize', 8,...
    'edgecolor', 'none', 'Position', [0 0.86 1 0.07], 'HorizontalAlignment', 'center' );

% add histogram
if histPlot == 1
    sb2=subplot(1,2,2);
    upL  = ceil(prctile(data(:),95));
    lowL = floor(prctile(data(:),5));
    data(data>upL)  = upL;
    data(data<lowL) = lowL;
    histogram(data,100, 'EdgeColor', 'none', 'FaceColor', [0 0 0])
    set(sb2, 'Box', 'off', 'Color', 'none', 'XLim', [lowL upL], 'Fontsize', 7);
    sb2.YAxis.Visible = 'off';
    
    sb1.Position = [0 0.12 1 0.75];
    sb2.Position = [0.14 0.33 0.15 0.17];
end

% savename
sp=['Map_' char(name) '.png'];

end


