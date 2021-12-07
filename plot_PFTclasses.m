%%%% plots the PFT classes
%% set save pth
spth = [pwd '/plots/PFT/'];
if ~exist(spth, 'dir'), mkdir(spth),end

%% load the PFT classification
f_in = 'data/input/globalTWS_Forcing.mat';
load(f_in, 'pft');

%% space stuff
load('data/input/lat_lon_global.mat', 'lat', 'lon')

%% map it
colTmp = othercolor('Paired12', 12);
colLabel.cticks  = [0:11];
colLabel.clabels = {'Sea', 'Ice', 'BEF', 'BDF', 'MixedF', 'CF', 'DF', 'WGrass', 'Shrubs', 'Tundra', 'Cult', 'Desert'}; %0:11; 
figure, PlotMapGlobal2('PFT classes ',[],[],lat, lon, pft, [-0.5 11.5],colTmp,colLabel); 
print(gcf,[spth 'PFT_classes.png'],'-dpng','-r300');

