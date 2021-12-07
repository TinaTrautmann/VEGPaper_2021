%%%%%%%% plot Koeppen Geiger zones
%% set save pth
spth = [pwd '/plots/StudyArea/'];
if ~exist(spth, 'dir'), mkdir(spth),end

% load space stuff
load('data/input/lat_lon_global.mat', 'lat', 'lon')

% load Koeppen Geiger zones
load('data/input/KoeppenGeiger_1deg.mat', 'KG');
[pix_x,pix_y] = LatLon2PixelCoord(lon,lat,90,-180,1);
idx = sub2ind([180, 360],pix_y, pix_x);

KG_map   = KG.Agg.KG_new';

figure, imagesc(KG_map)
KG_v          = KG_map(idx); %valids in map
KG_v(KG_v==1) = 2;
KG_v          = KG_v-1;

zonesID     = unique(KG_v);
zoneNames   = ['Trop' KG.Agg.CodesAgg(3:end)];

colLimKG    = [0.5 7.5];
colKG       = [rgb('DarkGreen');rgb('Gold'); rgb('YellowGreen');rgb('Olive');rgb('DarkCyan');rgb('CadetBlue');rgb('MediumBlue')];

colLabelKG.cticks  = [1:1:7];
colLabelKG.clabels  = zoneNames;

figure, PlotMapGlobal_noInfo('Koeppen Geiger Zones',[],[],lat,lon,KG_v,colLimKG,colKG,colLabelKG,1,1);
print(gcf,[spth 'Koeppen_Geiger_zones.png'],'-dpng','-r300');
