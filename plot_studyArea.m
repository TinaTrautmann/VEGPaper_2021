%%%% plot the study area and removed grids classes
%% set save pth
spth = [pwd '/plots/StudyArea/'];
if ~exist(spth, 'dir'), mkdir(spth),end

%% load the PFT classification
f_in = 'data/input/globalData.mat';
load(f_in, '');

%% global data
load(f_in, 'water_frac', 'permsnow_frac', 'artificial_frac', 'bare_frac', 'dirHuman', 'lat', 'lon')

% get the index in a 1deg map
[pix_x,pix_y] = LatLon2PixelCoord(lon,lat,90,-180,1);
idx = sub2ind([180, 360],pix_y, pix_x);

% checks
idxValid = find(water_frac<50 & permsnow_frac<10 & artificial_frac<10 & bare_frac<20 & isnan(dirHuman));

pix_v = length(idxValid);
lat_v = lat(idxValid);
lon_v = lon(idxValid);

%% show which grids have been extracted why
map_why                        = ones(size(lat));
map_why(water_frac>=50)        = 2;
map_why(permsnow_frac>=10)     = 3;
map_why(artificial_frac>=10)   = 4;
map_why(bare_frac>=20)         = 5;
map_why(~isnan(dirHuman))      = 6;

colTmp = [rgb('ForestGreen');...
            rgb('MediumBlue');...
            rgb('Turquoise');...
            rgb('Black');...
            rgb('SaddleBrown');...
            rgb('Purple')];
colLabel.cticks  = [1:6];
colLabel.clabels = {'study area', 'water', 'permsnow', 'artificial', 'bare', 'human impact TWS'};
figure, PlotMapGlobal2('Removed Grids from Study Area',[],[],lat, lon, map_why, [0.5 6.5],colTmp,colLabel); %'2 = water>50%, 3 = permsnow>10%, 4 = artificial>10%, 5 = bare>20%, 6 = direct human TWS trend'

%% load the cluster regions
load('data/input/clusterRegions.mat', 'CLregions');

KG_map   = CLregions;

KG_v        = KG_map(idx); %valids in map
KG_v(KG_v==4) = 10; 
KG_v(KG_v==2) = 20; 
KG_v(KG_v==3) = 30; 
KG_v(KG_v==1) = 40; 
KG_v(KG_v==5) = 50; 
KG_v(KG_v==0) = 4; 
KG_v = KG_v ./ 10;
KG_v(KG_v==0.4)=NaN;

zonesID     = unique(KG_v(~isnan(KG_v)));
zoneNames   = {'R1- Cold', 'R2- Temperate', 'R3- Humid',  'R4- Sub-humid',   'R5- Semi-arid'};
colLabel.clabels = zoneNames;
colLabel.cticks  = [1:5];
colKG = [rgb('DarkCyan');rgb('YellowGreen');rgb('DarkGreen');rgb('Olive');rgb('Gold')];
figure, PlotMapGlobal2('Cluster regions',[],[],lat, lon, KG_v, [0.5 5.5],colKG,colLabel); 
print(gcf,[spth 'StudyArea_ClusterRegions.png'],'-dpng','-r300');

%% combine the 2 maps
map_all = map_why+4;
map_all(map_why==1) = KG_v(map_why==1);

unique(map_all(~isnan(map_all)))
colLabel.clabels = {'R1', 'R2', 'R3', 'R4', 'R5', ...
                'w', 's', 'a', 'b', 'hTWS'};
colKG       = [rgb('DarkCyan');rgb('ForestGreen');rgb('DarkGreen');rgb('Olive');rgb('GoldenRod');...
                rgb('LightBlue'); rgb('White'); rgb('DarkGray'); rgb('NavajoWhite'); rgb('Thistle')];

colLabel.cticks  = [1:10];

figure, PlotMapGlobal2('Cluster regions of the study area & removed grids',[],[],lat, lon, map_all, [0.5 10.5],colKG, colLabel); 
print(gcf,[spth 'StudyArea_ClusterRegions_removedGrids.png'],'-dpng','-r300');


%% Percent Area
%% study area vs global 
study = load('data/input/lat_lon_global.mat');
glob  = load('data/input/globalData.mat', 'lat', 'lon');

study.area = AreaGridLatLon(study.lat, study.lon, [1 1]);
glob.area  = AreaGridLatLon(glob.lat, glob.lon, [1 1]);

study.area = study.area(:,1);
glob.area  = glob.area(:,1);

% study area vs global 
G = sum(glob.area(:));
P = sum(study.area(:));

perc = P/G;
disp(['study area vs. global land area: ' num2str(perc) ])

%% opti904
cal = load('data\input\lat_lon_opti904.mat');

cal.area = AreaGridLatLon(cal.lat, cal.lon, [1 1]);
cal.area = cal.area(:,1);

% opti vs global 
G = sum(glob.area(:));
P = sum(cal.area(:));

perc = P/G;
disp(['opti area vs. global land area: ' num2str(perc) ])


% opti vs study area
G = sum(study.area(:));
P = sum(cal.area(:));

perc = P/G;
disp(['opti area vs. study area: ' num2str(perc) ])

