function [ sp ] = PlotMapGlobal2(name, addinfo, units, lat, lon, data, colLim, col, colLabel)
% plots data onto a Map of the Northern Hemisphere, showing also land
% areas, large rivers and the study area
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


%
land    = shaperead('landareas', 'UseGeoCoords', true);
rivers  = shaperead('worldrivers', 'UseGeoCoords', true);
% area    = shaperead('C:\Profiles\ttraut\Dropbox\PLOTS\GeoData\valid_optims_shape\valid_optim_new.shp','UseGeoCoords',true);

geoRaRef    = georasterref('RasterSize', [180 360], 'RasterInterpretation', 'cells',  ...
            'LatitudeLimits', [-90 90], 'LongitudeLimits', [-180 180]);
Z           = NaN(180, 360); % prepare mapgrid
Z           = imbedm(lat, lon, data, Z, geoRaRef); % insert data in mapgrid
Z           = [Z ones(180,1)];  % add column for correct plotting of last column
Z           = [Z; ones(1,361)]; % add row for correct plotting of last row

if isempty(colLim)==1
    colLim  = [floor(min(data(:))) ceil(prctile(data(:),95))];
end

if isempty(col)==1
    col     = jet;
end

set(gcf,'Position', [0.5 0.5 16 9])
ax=worldmap('World'); 
set(ax, 'layer', 'top','xcolor', 'w', 'ycolor', 'w', 'Position', [0.1 0.12 0.8 0.75]);
axis off; grid on; axis tight;
setm(gca,'ParallelLabel', 'off','MeridianLabel','off');
geoshow(ax, land, 'FaceColor', [0.85 0.85 0.85])

if isfield(colLabel,'size')
    geoshow(ax, rivers, 'Color', [.25 .25 .25]);
    sm = scatterm(lat, lon, colLabel.size, data, 'filled')
else
    sm = surfm([-90 90], [-180 180], Z);
    geoshow(ax, rivers, 'Color', [.25 .25 .25])
end

hold on;
% geoshow(ax, area, 'EdgeColor', [.7 0 0], 'LineWidth', 0.5, 'FaceColor', [1 1 1], 'FaceAlpha', 0.1)

try
    caxis(colLim)
end;
colormap(col)
cb  =   colorbar;
cb1 =   ylabel(cb, units);
t   =   title(strrep(name,'_','-'));

% edit colorbar
set(cb,'FontSize',8, 'Position', [0.1 0.09 0.8 0.02],...
    'Orientation', 'horizontal') 
try  % set limits and ticknames if they are specified
    cticks  = colLabel.cticks;
    clabels = colLabel.clabels;
    set(cb, 'YTick', cticks, 'YTickLabel', clabels);
end
set(cb1,'Position',[(colLim(2)+colLim(1))/2 -2.75 0])

% title and heading
set(t,'Fontsize', 9, 'Fontweight', 'b', 'Position', [9 10000000 0 ]); %worldmap
%set(t,'Fontsize', 9, 'Fontweight', 'b', 'Position', [9 10000000 0 ]); %EU map
annotation('textbox', 'String', addinfo, 'fontsize', 8,...
    'edgecolor', 'none', 'Position', [0 0.85 1 0.07], 'HorizontalAlignment', 'center' );

pos= get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Centimeter','PaperSize',[pos(3), pos(4)])

% savename
sp=['Map_' char(name)];


end

