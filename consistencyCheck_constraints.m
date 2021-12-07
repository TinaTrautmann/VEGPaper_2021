%% ------------ WB check of constraints for VEG Paper ------------ %%
%% set save pth
spth = [pwd '/plots/consistency_check_obs/'];
if ~exist(spth, 'dir'), mkdir(spth),end


%% load the obs
obs_in = 'data/input/globalBaseline_Constraints_1deg.mat';
f_in   = 'data/input/globalTWS_Forcing.mat';

% observational constraints
load(obs_in, 'TWS_GRACE');
TWS   = TWS_GRACE; clear TWS_GRACE;

load(obs_in, 'ET_Fluxcom');
Evap  = ET_Fluxcom; clear ET_Fluxcom;

load(obs_in, 'Q_GRUN');
Qrun  = Q_GRUN; clear Q_GRUN;

% forcing
load(f_in, 'Rain');
load(f_in, 'Snow');

% coordinates
load('data/input/lat_lon_global.mat', 'lat', 'lon')

% scaling parameter for snowfall
p_SF_scale    = 0.67; % from Trautmann et al. 2018

% time
time_data = {'2000-03-01', '2014-12-31'}; % of model run
time_WB   = {'2004-01-01', '2010-11-30'}; % for consistency check, without any NaNs

[idx_data, idx_WB, time, xMonth] = calcConsisTime(time_data, time_WB, 'monthly');


%% Data
% daily values to monthly means
Precip_d = Rain + p_SF_scale .* Snow; 
Precip   = aggDay2Mon(Precip_d, time_data{1}, time_data{2});  

% adjust time
idx_data_TWS = [idx_data idx_data(end)+1];

P       = Precip(:,idx_data).*30;
ET      = Evap(:,idx_data).*30;
Q       = Qrun(:,idx_data).*30;

TWS_0   = TWS(:,idx_data_TWS);
m       = mean(TWS_0,2);
m       = repmat(m,1,size(TWS_0,2)); %recalculate time-mean anomaly

TWS_0   = TWS_0-m+10000; %add 10000mm to make it all positive

% monthly storage changes
dS = diff(TWS_0,[],2);


%% calculate monthly WB
WB0      = P - ET - Q - dS;

P_WB0    = ET + Q + dS;
ET_WB0   = P - Q - dS;
Q_WB0    = P - ET - dS;
dS_WB0   = P - ET - Q;


%% Figures
% Map WB
figure
sp = PlotMapGlobal2('Mean Water Imbalance',[time{1} ' to ' time{2}],'mm',lat, lon, nanmean(WB0,2),[-100 100],othercolor('BuDRd_18',21),[])
print(gcf, [spth 'meanWB_100-100.png'], '-dpng', '-r300')

figure
sp = PlotMapGlobal2('Median Water Imbalance',[time{1} ' to ' time{2}],'mm',lat, lon, nanmedian(WB0,2),[-100 100],othercolor('BuDRd_18',41),[])
print(gcf, [spth 'medianWB_100-100.png'], '-dpng', '-r300')


% WB scaled by Precip
P_mean        = nanmean(P,2);
WB_mean       = nanmean(WB0,2); 
WB_scaled     = WB_mean ./ P_mean;
 

figure
sp = PlotMapGlobal_noInfo(['Mean Water Imbalance (scaled by P)'],['mean = ' num2str(round(nanmean(WB_scaled(:)),2)) ' | median = ' num2str(round(median(WB_scaled(:),'omitnan'),2))],'',lat, lon, WB_scaled, [-1 1],othercolor('BuDRd_18',250),[],1,1)
print(gcf, [spth 'meanWB_scaled_statistics.png'], '-dpng', '-r300')


%% WB time series in in 1 plot
WBcomps = {'P', 'ET', 'Q', 'dS'};

% time series with same y-limits
yLim    = [0 90] 
ylim_t  = yLim(2)-yLim(1);
ylim_t  = ylim_t /2;
ylim_t  = [0-ylim_t ylim_t];
    
figure, set(gcf, 'Position', [5 2 16 20])
ha = tight_subplot(5,1,[.05 .05],[.03 .07],[.07 .05]);
axes(ha(1))            
plot(xMonth, nanmean(WB0,1), 'b-')
hold on, plot([xMonth(1) xMonth(end)],[0 0], '-k')
title(['Monthly Water Imbalance from observations ' ]);
ylabel('mm'), set(gca, 'fontsize', 8,'YLim', ylim_t);
legend('P - ET - Q - dS', 'Location', 'SouthOutside', 'Orientation', 'Horizontal');
for wb = 1: numel(WBcomps)   
    compName = WBcomps{wb};

    dataObs  = eval(char(compName));

    dataWB = eval(char([compName '_WB0']));
    tmpObs = nanmean(dataObs,1);
    tmpWB  = nanmean(dataWB,1);
           
    axes(ha(wb+1))
    plot(xMonth, tmpObs, 'r-')
    hold on, plot(xMonth, tmpWB, 'b-'),
    title(['Monthly ' compName ])
    legend({'obs', ['WB- r = ' num2str(round(corr(tmpObs(:), tmpWB(:)),2)) ]}, 'Location', 'SouthOutside', 'Orientation', 'Horizontal')
    ylabel('mm'), set(gca, 'fontsize', 8)
    if strcmp(compName,'dS') 
        set(gca,'YLim', ylim_t);
    else
        set(gca,'YLim', yLim);
    end
end
hold on, plot([xMonth(1) xMonth(end)],[0 0], '-k')
print(gcf, [spth 'WB_from_Constraints_sameYaxis.png'], '-dpng', '-r300')


%% Load & add the modelled data of B & VEG
expNames    = {'B', 'E_B_bL_RD4'}; 
expNames2   = {'B', 'VEG'};
rDates      = {'20200629', '20200629'};

varNames = {'wTotal', 'evapTotal', 'roTotal', 'evapSub'};

%loop over experiments
for eN =1:numel(expNames)
    expName = expNames{eN};
    rDate   = rDates{eN};
    % path of model output
    exp_in  = ['data/model_output/validCali_' expName '_validCali_' rDate '/modelOutput/'];
    
    for vn=1:numel(varNames)
        vName = varNames{vn};
        fName = [exp_in 'validCali_' expName '_validCali_' vName '.nc'];
        if isfile(fName)
            mod.(expName).(vName)          =    ncread(fName,vName);
        end
    end
    
    %load info
    infoIn.(expName) = load([exp_in 'validCali_' expName '_' rDate '_info.mat']);
end

% process and adjust the data
% use the monthly sum of Precip, instead of average *30
P_mod = sumDay2Mon(Precip_d, time_data{1}, time_data{2});
P_mod = P_mod(:,idx_data);

for eN =1:numel(expNames)
    expName = expNames{eN};
    for vn=1:numel(varNames)
        vName = varNames{vn};
        % extract the consistent time
        if strcmp(vName,'wTotal')
            mod.(expName).(vName) = squeeze(mod.(expName).(vName));
            % make monthly means
            tmp_monthly = aggDay2Mon(mod.(expName).(vName), time_data{1}, time_data{2});
            
            % time mean
            tmp_monthly_ano   = tmp_monthly(:,idx_data_TWS);
            m       = mean(tmp_monthly_ano,2);
            m       = repmat(m,1,size(tmp_monthly_ano,2));            
            tmp_monthly_ano   = tmp_monthly_ano-m+10000;
            mod_month.(expName).(vName)     = tmp_monthly_ano;
            % monthly storage changes
            mod_month.(expName).dS          = diff(mod_month.(expName).(vName),[],2);
        elseif strcmp(vName,'evapTotal')
            % make monthly sums
            tmp_monthly = sumDay2Mon(mod.(expName).(vName), time_data{1}, time_data{2});            
            mod_month.(expName).ET          = tmp_monthly(:,idx_data);
        elseif strcmp(vName,'roTotal')
            % make monthly sums
            tmp_monthly = sumDay2Mon(mod.(expName).(vName), time_data{1}, time_data{2});
            mod_month.(expName).Q           = tmp_monthly(:,idx_data);
        elseif strcmp(vName,'evapSub')
            % make monthly sums
            tmp_monthly = sumDay2Mon(mod.(expName).(vName), time_data{1}, time_data{2});
            mod_month.(expName).ET_Sub      = tmp_monthly(:,idx_data);
       end
    end
    
    % sublimation is not included in the evapTotal -> but needs to be
    % considered for the WB
%     WB_mod.(expName) = P_mod - mod_month.(expName).ET - mod_month.(expName).Q - mod_month.(expName).ET_Sub - mod_month.(expName).dS;
    
    % do the test by calculating the daily WB
    tmp_WB           = calcWaterBalance(Precip_d,squeeze(mod.(expName).wTotal),mod.(expName).evapTotal+mod.(expName).evapSub+mod.(expName).roTotal);
    WB_mod.(expName) = sumDay2Mon(tmp_WB, '2000-03-02', time_data{2});
    WB_mod.(expName) = WB_mod.(expName)(:,idx_data);
end


%% add the mod to the time series
colExp  = [rgb('MediumBlue');rgb('ForestGreen')];
    
WBcomps = {'P', 'ET', 'Q', 'dS'};

figure, set(gcf, 'Position', [5 1 16 21])
ha = tight_subplot(5,1,[.05 .05],[.03 .07],[.07 .05])
axes(ha(1))            
plot(xMonth, nanmean(WB0,1), 'k-'), hold on
for eN =1:numel(expNames)
    plot(xMonth, nanmean(WB_mod.(expNames{1}),1), '-', 'color', colExp(eN,:));
end
hold on, plot([xMonth(1) xMonth(end)],[0 0], ':k')
title(['Monthly Water Imbalance | P - ET - Q - dS' ])
ylabel('mm'), set(gca, 'fontsize', 8)
tmpRMSE = calcRMSE(zeros(1,length(nanmean(WB0,1))), nanmean(WB0,1));
l = legend([['from obs | RMSE = ' num2str(round(tmpRMSE,1)) 'mm'], expNames2], 'Location', 'SouthOutside', 'Orientation', 'Horizontal');
l.Box = 'off';
tmpPos = l.Position;
tmpPosAx = ha(1).Position;
l.Position = [0.01 l.Position(2) 0.99 l.Position(4)]
ha(1).Position = tmpPosAx;
for wb = 1: numel(WBcomps)   
    compName = WBcomps{wb};
    dataObs  = eval(char(compName));

    dataWB = eval(char([compName '_WB0']));
    tmpObs = nanmean(dataObs,1);
    tmpWB  = nanmean(dataWB,1);
            
    axes(ha(wb+1));
    plot(xMonth, tmpObs, 'r:', 'LineWidth', 1.25)
    hold on, plot(xMonth, tmpWB, 'k-'),
    tmpRMSE = calcRMSE(tmpObs(:), tmpWB(:));
    leg      = {['from obs | r=' num2str(round(corr(tmpObs(:), tmpWB(:)),2)) ' ; RMSE=' num2str(round(tmpRMSE,1)) 'mm']};
    if strcmp(compName,'P')
%         plot(xMonth, nanmean(P_mod,1), '-', 'color', colExp(1,:))
    else
        for eN =1:numel(expNames)
            tmpExp = nanmean(mod_month.(expNames{eN}).(compName),1);
            plot(xMonth,tmpExp , '-', 'color', colExp(eN,:))
            tmpRMSE = calcRMSE(tmpObs(:), tmpExp(:));
            leg = [leg [expNames2{eN} ' | r=' num2str(round(corr(tmpObs(:), tmpExp(:)),2)) ' ; RMSE=' num2str(round(tmpRMSE,1)) 'mm']];
        end
    end
    title(['Monthly ' compName ])
    ylabel('mm'), set(gca, 'fontsize', 8)
    l=legend(['obs', leg], 'Location', 'SouthOutside', 'Orientation', 'Horizontal', 'fontsize',7);
    l.Box = 'off';
    
    tmpPos = l.Position;
    tmpPosAx = ha(wb+1).Position;
    l.Position = [0.01 l.Position(2) 0.99 l.Position(4)]
    ha(wb+1).Position = tmpPosAx;
end
% hold on, plot([xMonth(1) xMonth(end)],[0 0], '-k')
print(gcf, [spth 'WB_from_Constraints_vs_model_rmse.png'], '-dpng', '-r300')

%% MSC regions - model-obs vs WB-obs
[~, ~, ~, ~, ~, M] = createDateVector(time{1}, time{2}, 'm');

pix_a = AreaGridLatLon(lat,lon,[1 1]);
pix_a = pix_a(:,1);

% get the zones
load('data/input/clusterRegions.mat', 'CLregions');
[pix_x,pix_y] = LatLon2PixelCoord(lon,lat,90,-180,1);
idx = sub2ind([180, 360],pix_y, pix_x);

KG_map   = CLregions;

KG_v        = KG_map(idx); %valids in map
KG_v(KG_v==0) = 4; 
KG_v(KG_v==4) = 10; 
KG_v(KG_v==2) = 20; 
KG_v(KG_v==3) = 30; 
KG_v(KG_v==1) = 40; 
KG_v(KG_v==5) = 50; 
KG_v = KG_v ./ 10;

zonesID     = unique(KG_v);
zoneNames   = {'R1- Cold', 'R2- Temperate', 'R3- Humid',  'R4- Sub-humid',   'R5- Semi-arid'};


% % use Koeppen Geiger zones
% load('data/input/KoeppenGeiger_1deg.mat', 'KG');
% [pix_x,pix_y] = LatLon2PixelCoord(lon,lat,90,-180,1);
% idx = sub2ind([180, 360],pix_y, pix_x);
% 
% KG_map   = KG.Agg.KG_new';
% 
% KG_v        = KG_map(idx); %valids in map
% KG_v(KG_v==1) = 2;
% KG_v        = KG_v-1;
% 
% zonesID     = unique(KG_v);
% zoneNames   = ['Trop' KG.Agg.CodesAgg(3:end)];
% 

unitNames = {'mm', 'mm', 'mm'};

WBcomps2  = WBcomps(2:end);
for wb = 1: numel(WBcomps2)   
    compName = WBcomps2{wb};
    
    %obs
    dataObs     = eval(char(compName));
    oData_MSC   = calcMSC(dataObs, M);
    
    %WB
    wbData      = eval(char([compName '_WB0']));
    wbData_MSC  = calcMSC(wbData, M);
    WBdiffObsMSC.(compName)   = wbData_MSC- oData_MSC;
    
    %mod
    for eN =1:numel(expNames)
        expName = expNames{eN};
        sData       = mod_month.(expName).(compName);
        sData_MSC   = calcMSC(sData, M);
        eval(char(['mod' num2str(eN) 'diffObsMSC.(compName)   = sData_MSC - oData_MSC;']))
    end
end
    

%% Plot the MSC for zones -with RMSE to zero
[sname, T_rmse0, sp] = plotMSCvarsZones_rmse0({'WB - obs', 'B - obs', 'VEG - obs'}, WBdiffObsMSC, mod1diffObsMSC, mod2diffObsMSC, zoneNames, KG_v, pix_a, unitNames, []);
print(gcf,[spth sname '_sameYLim_RMSE.png'],'-dpng','-r300');


