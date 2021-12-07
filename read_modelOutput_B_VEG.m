
%% set save pth
spth = [pwd '/plots/B_vs_VEG/'];
if ~exist(spth, 'dir'), mkdir(spth),end


%% load model output
expNames    = {'B', 'E_B_bL_RD4'}; 
rDates      = {'20200629', '20200629'};

space_run = 'validCali';

varNames = {'wTotal' 'evapTotal', 'roTotal', 'wSoil', 'wGW', 'wSurf', 'wSnow',...
                'tranAct','evapSoil', 'evapInt', 'roSurfIndir', 'roSurfDir'};

            
%loop over experiments
for eN =1:numel(expNames)
    expName = expNames{eN};
    rDate   = rDates{eN};
    % path of model output
    exp_in  = ['data/model_output/validCali_' expName '_validCali_' rDate '/modelOutput/'];
    % load variables
    for vn=1:numel(varNames)
        vName = varNames{vn};
        fName = [exp_in 'validCali_' expName '_validCali_' vName '.nc'];
        if isfile(fName)
            mod.(expName).(vName) = ncread(fName,vName);
        end
    end
    
    %load info
    infoIn.(expName) = load([exp_in 'validCali_' expName '_' rDate '_info.mat']);
end

info = infoIn.(expNames{1}).info;

%% load observations 
obs_in = 'data/input/globalBaseline_Constraints_1deg.mat'

obsNames_in = {'TWS_GRACE', 'wSnow', 'ET_Fluxcom', 'wSoil', 'Q_GRUN'};
uncNames_in = {'TWS_GRACE_unc', 'wSnow_unc', 'ET_Fluxcom_unc', 'wSoil_unc', 'Q_GRUN_unc'}; 
obsNames_an = {'TWS', 'SWE', 'Evap', 'wSoil', 'Q'};

for oi=1:numel(obsNames_an)
   tmpData  = load(obs_in, obsNames_in{oi});
   tmpUnc   = load(obs_in, uncNames_in{oi});
   obs.(obsNames_an{oi}).data = tmpData.(obsNames_in{oi});
   obs.(obsNames_an{oi}).unc  = tmpUnc.(uncNames_in{oi});
end

% %% load forcing
% f_in = 'data/input/globalTWS_Forcing.mat'
% f    = load(f_in, 'Rn', 'Rain', 'Snow', 'Tair', 'PSurfDay', 'PET') %thats just the basic, more data is needed for the VEG experiment
 
%% -- Prep Stuff --
% space stuff
load('data/input/lat_lon_global.mat', 'lat', 'lon')
pix_a = AreaGridLatLon(lat,lon,[1 1]);
pix_a = pix_a(:,1);

% time stuff
[xDay,  ~, ~, ~, ~, Md,D]   = createDateVector(info.tem.model.time.sDate,info.tem.model.time.eDate, 'd');
[xMonth, ~, ~, ~, ~, M]     = createDateVector(info.tem.model.time.sDate,info.tem.model.time.eDate, 'm');

tick_locations  = xMonth(M==1);
tick_locationsD = xDay(D==1);

nPix   = size(pix_a,1);
nTix   = size(xDay,2);

% colors
colLabel = [];


%% get the zones
% load the cluster regions
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

colLimKG    = [0.5 5.5];
colKG       = [rgb('DarkCyan');rgb('YellowGreen');rgb('DarkGreen');rgb('Olive');rgb('Gold')];
colLabelKG  = colLabel;
colLabelKG.cticks  = [1:1:5];
colLabelKG.clabels  = zoneNames;

figure, PlotMapGlobal_noInfo('Cluster derived Regions',[],[],lat,lon,KG_v,colLimKG,colKG,colLabelKG,1,1);
%      print([spth 'clusterRegions.png'],'-dpng','-r300')

% % use Koeppen Geiger zones
% load('data/input/KoeppenGeiger_1deg.mat', 'KG');
% [pix_x,pix_y] = LatLon2PixelCoord(lon,lat,90,-180,1);
% idx = sub2ind([180, 360],pix_y, pix_x);
% 
% KG_map   = KG.Agg.KG_new';
% 
% figure, imagesc(KG_map)
% KG_v        = KG_map(idx); %valids in map
% KG_v(KG_v==1) = 2;
% KG_v        = KG_v-1;
% 
% zonesID     = unique(KG_v);
% zoneNames   = ['Trop' KG.Agg.CodesAgg(3:end)];
% 
% 
% colLimKG    = [0.5 7.5];
% colKG       = [rgb('DarkGreen');rgb('Gold'); rgb('YellowGreen');rgb('Olive');rgb('DarkCyan');rgb('CadetBlue');rgb('MediumBlue')];
% 
% colLabelKG  = colLabel;
% colLabelKG.cticks  = [1:1:7];
% colLabelKG.clabels  = zoneNames;
% 
% figure, PlotMapGlobal_noInfo('Koeppen Geiger Zones',[],[],lat,lon,KG_v,colLimKG,colKG,colLabelKG,1,1);
% print(gcf,[spth 'Koeppen_Geiger_zones.png'],'-dpng','-r300');

%


%% --- FIGURES
% only consider TWSobs >= -500mm and <= 500mm
TWScat = 1; 

optNames = expNames;

% savepaths
spth_p      = [spth 'ModelPerformance/' ];
if ~exist(spth_p, 'dir'), mkdir(spth_p),end

% what variables to compare
dataNames    = ['obs', optNames];
names       = {'wTWS', 'ET', 'Q'} ;
sDatas      = {'reshape(mod.XX.wTotal, info.tem.helpers.sizes.nPix, info.tem.helpers.sizes.nTix)',...
    'mod.XX.evapTotal', 'mod.XX.roTotal'};
oDatas      = {'obs.TWS.data', 'obs.Evap.data', 'obs.Q.data'};
oDatas_unc  = {'obs.TWS.unc',  'obs.Evap.unc', 'obs.Q.unc'};
unitNames   = {'mm', 'mm/d', 'mm/d'};

% Loop over variables & put everything into structures for obs, exp1, exp2
for ii=1:numel(names)
    name  = names{ii};
    oData = eval([oDatas{ii}]);
    
    % anomaly if TWS
    if strcmp(names{ii},'wTWS')
        if TWScat == 1
            oData(oData<=-500) = NaN;
            oData(oData>=500)  = NaN;
        end
        m = nanmean(oData,2);
        oData = oData - repmat(m,1,length(xMonth));
    end
    NaNidx  = find(isnan(oData));
        
    % MSC
    oData_MSC   = calcMSC(oData, M);
    
    % put into structure
    dataObs.(name)      = oData;
    dataObsMSC.(name)   = oData_MSC;
    
    % loop over experiments
    for op=1:numel(optNames)
        optN = optNames{op};
        tmp = strrep(sDatas{ii}, 'XX', [optN]);
        sData_d = eval(tmp);
        % monthly aggregation
        sData = aggDay2Mon(sData_d, info.tem.model.time.sDate, info.tem.model.time.eDate);
        % only consistent data points
        sData(NaNidx) = NaN;
        
        % anomaly if TWS
        if strcmp(names{ii},'wTWS')
            m = nanmean(sData,2);
            sData = sData - repmat(m,1,length(xMonth));
        end
        % calculate MSC
        sData_MSC   = calcMSC(sData, M);
        
        % put into structure
        eval(char(['data' num2str(op) '.(name)      = sData;']))
        eval(char(['data' num2str(op) 'MSC.(name)   = sData_MSC;']))
    end
end

%% Plot the MSC for zones
[sname, T_corr, sp] = plotMSCvarsZones_VEG1(dataNames, dataObsMSC, data1MSC, data2MSC, zoneNames, KG_v, pix_a, unitNames, []);
print(gcf,[spth_p sname '.png'],'-dpng','-r300');


%% Plot the Correlation Maps
plotDiffMaps2(dataNames, dataObsMSC, data1MSC, data2MSC, lat, lon, [], spth_p)
close all


%% --- calculate global metrics
dataUnc.wTWS    = obs.TWS.unc;
dataUnc.ET      = obs.Evap.unc;
dataUnc.Q       = obs.Q.unc;
dataObs.wTWS    = obs.TWS.data;
dataObs.ET      = obs.Evap.data;
dataObs.Q       = obs.Q.data;

% loop over variables
varNames    = {'wTWS','ET','Q'}
varNamesExp = {'wTotal','evapTotal','roTotal'};
for vN=1:numel(varNames)
    varN  = varNames{vN};
    varNE = varNamesExp{vN};
    
    T_V = array2table(NaN(numel(expNames),6));
    T_V.Properties.RowNames        = expNames;

    % loop over experiments
    for eN=1:numel(expNames)
        expN = expNames{eN};
        
        % get the data
       if ndims(mod.(expN).(varNE))>2
           tmp_data =  sum(mod.(expN).(varNE),2);
           data.(expN).(varN) = reshape(tmp_data, info.tem.helpers.sizes.nPix, info.tem.helpers.sizes.nTix);
       else
           data.(expN).(varN) = mod.(expN).(varNE);
       end
        
       % monthly aggregation
       data.(expN).(varN) = aggDay2Mon(data.(expN).(varN), info.tem.model.time.sDate, info.tem.model.time.eDate);
       
        % anomaly if TWS
        if strcmp(varN,'wTWS')
            if TWScat == 1
                dataObs.(varN)(dataObs.(varN)<=-500) = NaN;
                dataObs.(varN)(dataObs.(varN)>=500)  = NaN;
            end
            m = nanmean(dataObs.(varN),2);
            dataObs.(varN) = dataObs.(varN) - repmat(m,1,length(xMonth));
            
            NaNidx  = find(isnan(dataObs.(varN)));
            data.(expN).(varN)(NaNidx) = NaN;
            m = nanmean(data.(expN).(varN),2);
            data.(expN).(varN) = data.(expN).(varN) - repmat(m,1,length(xMonth));
        else
            NaNidx  = find(isnan(dataObs.(varN)));
            data.(expN).(varN)(NaNidx) = NaN;
        end

        % calculate metrics
        % KGE (only calculated for ~isnan)
        [KGE, r_p, alpha, beta] = calcKGE(dataObs.(varN)(:),data.(expN).(varN)(:));
        S_VE.KGE      = KGE;
        S_VE.corr     = r_p;
        S_VE.var_err  = (alpha-1)^2;
        S_VE.bias_err = (beta-1)^2;        
        
        % MEF
        v_obs   = find(~isnan(dataObs.(varN)(:)));
        NSE     = calcMEF(dataObs.(varN)(v_obs),data.(expN).(varN)(v_obs), ones(size(data.(expN).(varN)(v_obs))));
        S_VE.MEF    = 1-NSE;
        
        % weighted MEF
        wNSE        = calcMEF(dataObs.(varN)(v_obs),data.(expN).(varN)(v_obs), dataUnc.(varN)(v_obs));
        S_VE.wMEF   = 1-wNSE;
         
        tmp_T         = struct2table(S_VE);
        T_V(expN,:)   = tmp_T;
        
        % metric names
        T_V.Properties.VariableNames = fieldnames(S_VE);
        
    end
        % output table
        
        eval(char(['T_' varN '= T_V;']))
end

% write the output
writetable(T_wTWS,[spth space_run '_globalPerformance.xls'],'Sheet','wTWS','WriteRowNames',true)
writetable(T_ET,[spth space_run '_globalPerformance.xls'],'Sheet','ET','WriteRowNames',true)
writetable(T_Q,[spth space_run '_globalPerformance.xls'],'Sheet','Q','WriteRowNames',true)

% same for MSC
for vN=1:numel(varNames)
    varN  = varNames{vN};
    varNE = varNamesExp{vN};
    
    T_V = array2table(NaN(numel(expNames),6));
    T_V.Properties.RowNames        = expNames;

    % loop over experiments
    for eN=1:numel(expNames)
        expN = expNames{eN};
        
        
       % MSC
        dataObsMSC.(varN)       = calcMSC(dataObs.(varN), M);
        dataMSC.(expN).(varN)   = calcMSC(data.(expN).(varN), M);
        dataUncMSC.(varN)       = calcMSC(dataUnc.(varN), M);
        % calculate metrics
        % KGE (only calculated for ~isnan)
        [KGE, r_p, alpha, beta] = calcKGE(dataObsMSC.(varN)(:),dataMSC.(expN).(varN)(:));
        S_VE.KGE      = KGE;
        S_VE.corr     = r_p;
        S_VE.var_err  = (alpha-1)^2;
        S_VE.bias_err = (beta-1)^2;        
        
        % MEF
        v_obs   = find(~isnan(dataObsMSC.(varN)(:)));
        NSE     = calcMEF(dataObsMSC.(varN)(v_obs),dataMSC.(expN).(varN)(v_obs), ones(size(dataMSC.(expN).(varN)(v_obs))));
        S_VE.MEF    = 1-NSE;
        
        % weighted MEF
        wNSE     = calcMEF(dataObsMSC.(varN)(v_obs),dataMSC.(expN).(varN)(v_obs), dataUncMSC.(varN)(v_obs));
        S_VE.wMEF    = 1-wNSE;
        
        tmp_T         = struct2table(S_VE);
        T_V(expN,:)   = tmp_T;
        
        % metric names
        T_V.Properties.VariableNames = fieldnames(S_VE);
        
    end
        % output table
        
        eval(char(['T_' varN '= T_V;']))
end

% write the output
writetable(T_wTWS,[spth space_run '_globalPerformance_MSC2.xls'],'Sheet','wTWS','WriteRowNames',true)
writetable(T_ET,[spth space_run '_globalPerformance_MSC2.xls'],'Sheet','ET','WriteRowNames',true)
writetable(T_Q,[spth space_run '_globalPerformance_MSC2.xls'],'Sheet','Q','WriteRowNames',true)


%% TWS composition as Impact Index (Getirana et al. 2017) - all storages seperate %%
withObs     = 0;
if withObs==1
    spth_p      = [spth 'TWScomp_validObs/' ];
else
    spth_p      = [spth 'TWScomp/' ];
end
if ~exist(spth_p, 'dir'), mkdir(spth_p),end

%prepare TWSobs
if withObs==1
    oData = obs.TWS.data;
    if TWScat == 1
        oData(oData<=-500) = NaN;
        oData(oData>=500)  = NaN;
    end
    idxNaN = find(isnan(oData));
    
    oData = DetrendMatrix(oData);
    %recalculate the TWS time mean
    m = nanmean(oData,2);
    oData = oData - repmat(m,1,length(xMonth));
    %calculate MSC & IAV
    [seasonObs,~,iavObs,~] = calcMSC(oData,M);
    %put into structure
    TS.msc_ano.obs.TWSobs = seasonObs;
    TS.iav_ano.obs.TWSobs = iavObs;
    TS.ano.obs.TWSobs     = oData;
end

% MSC - sum of absolute anomaly to MSC
varNames = {'IwSnow','IwSoil','IwGW', 'IwSurf'};
for eN=1:numel(expNames)
    expName = expNames{eN};
    
    wTotal_d = reshape(mod.(expName).wTotal, info.tem.helpers.sizes.nPix, info.tem.helpers.sizes.nTix);
    wSoil_d  = reshape(mod.(expName).wSoil(:,1,:), info.tem.helpers.sizes.nPix, info.tem.helpers.sizes.nTix) + reshape(mod.(expName).wSoil(:,2,:), info.tem.helpers.sizes.nPix, info.tem.helpers.sizes.nTix);
    wSnow_d  = reshape(mod.(expName).wSnow, info.tem.helpers.sizes.nPix, info.tem.helpers.sizes.nTix);
    wGW_d    = reshape(mod.(expName).wGW, info.tem.helpers.sizes.nPix, info.tem.helpers.sizes.nTix);
    wSurf_d  = reshape(mod.(expName).wSurf, info.tem.helpers.sizes.nPix, info.tem.helpers.sizes.nTix);

    % aggregate to months
    wTotal_m = aggDay2Mon(wTotal_d, info.tem.model.time.sDate, info.tem.model.time.eDate);
    wSnow_m  = aggDay2Mon(wSnow_d, info.tem.model.time.sDate, info.tem.model.time.eDate);
    wSoil_m  = aggDay2Mon(wSoil_d, info.tem.model.time.sDate, info.tem.model.time.eDate);
    wGW_m    = aggDay2Mon(wGW_d, info.tem.model.time.sDate, info.tem.model.time.eDate);
    wSurf_m  = aggDay2Mon(wSurf_d, info.tem.model.time.sDate, info.tem.model.time.eDate);
    
    % only validObs?
    if withObs==1
        wTotal_m(idxNaN) = NaN;
        wSnow_m(idxNaN) = NaN;
        wSoil_m(idxNaN) = NaN;
        wGW_m(idxNaN) = NaN;
        wSurf_m(idxNaN) = NaN;
    end

    % calculate MSC and IAV
    [wTotal_msc,~, wTotal_iav] = calcMSC(wTotal_m,M);
    [wSnow_msc,~, wSnow_iav]   = calcMSC(wSnow_m,M);
    [wSoil_msc,~, wSoil_iav]   = calcMSC(wSoil_m,M);
    [wGW_msc,~, wGW_iav]       = calcMSC(wGW_m,M);
    [wSurf_msc,~, wSurf_iav]   = calcMSC(wSurf_m,M);

    
    % ----for msc & iav
    aggTimes = {'msc', 'iav'};
    for ii=1:numel(aggTimes)
        aggTime = aggTimes{ii};
        
        wSnow  = eval(char(['wSnow_' aggTime ';']));
        wSoil  = eval(char(['wSoil_' aggTime ';']));
        wGW    = eval(char(['wGW_' aggTime ';']));
        wSurf  = eval(char(['wSurf_' aggTime ';']));
        wTWS   = eval(char(['wTotal_' aggTime ';']));
        
        % calculate variances & covariances - MSC
        C_wTotal   = NaN(nPix,1);
        C_wSnow    = NaN(nPix,1);
        C_wSoil    = NaN(nPix,1);
        C_wGW      = NaN(nPix,1);
        C_wSurf    = NaN(nPix,1);
        
        for i = 1:nPix
            C_wTotal(i,1)  = sum(abs(wTWS(i,:) - nanmean(wTWS(i,:))));
            C_wSnow(i,1)   = sum(abs(wSnow(i,:) - nanmean(wSnow(i,:))));
            C_wSoil(i,1)   = sum(abs(wSoil(i,:) - nanmean(wSoil(i,:))));
            C_wGW(i,1)     = sum(abs(wGW(i,:) - nanmean(wGW(i,:))));
            C_wSurf(i,1)   = sum(abs(wSurf(i,:) - nanmean(wSurf(i,:))));
        end
                
        CR.(expName).(aggTime).C_wSnow  = C_wSnow;        
        CR.(expName).(aggTime).C_wSoil  = C_wSoil;
        CR.(expName).(aggTime).C_wGW    = C_wGW;
        CR.(expName).(aggTime).C_wSurf  = C_wSurf;
        CR.(expName).(aggTime).C_wTotal = C_wTotal;
        
        % impact index
        CR.(expName).(aggTime).IwSnow  = C_wSnow./(C_wSnow+C_wSoil+C_wGW+C_wSurf);
        CR.(expName).(aggTime).IwSoil  = C_wSoil./(C_wSnow+C_wSoil+C_wGW+C_wSurf);
        CR.(expName).(aggTime).IwGW    = C_wGW./(C_wSnow+C_wSoil+C_wGW+C_wSurf);
        CR.(expName).(aggTime).IwSurf  = C_wSurf./(C_wSnow+C_wSoil+C_wGW+C_wSurf);
                
    end
    
    
       % ----Index of global average msc & iav
    aggTimes = {'msc', 'iav'};
    for ii=1:numel(aggTimes)
        aggTime = aggTimes{ii};
        
        wSnow  = eval(char(['nanmean(wSnow_' aggTime ',1);']));
        wSoil  = eval(char(['nanmean(wSoil_' aggTime ',1);']));
        wGW    = eval(char(['nanmean(wGW_' aggTime ',1);']));
        wSurf  = eval(char(['nanmean(wSurf_' aggTime ',1);']));
        wTWS   = eval(char(['nanmean(wTotal_' aggTime ',1);']));
        
        % calculate contribution
        C_wTotal  = sum(abs(wTWS - nanmean(wTWS(:))));
        C_wSnow   = sum(abs(wSnow - nanmean(wSnow(:))));
        C_wSoil   = sum(abs(wSoil - nanmean(wSoil(:))));
        C_wGW     = sum(abs(wGW - nanmean(wGW(:))));
        C_wSurf   = sum(abs(wSurf - nanmean(wSurf(:))));
        
        CR_Avg.(expName).(aggTime).C_wSnow  = C_wSnow;
        CR_Avg.(expName).(aggTime).C_wSoil  = C_wSoil;
        CR_Avg.(expName).(aggTime).C_wGW    = C_wGW;
        CR_Avg.(expName).(aggTime).C_wSurf  = C_wSurf;
        CR_Avg.(expName).(aggTime).C_wTotal = C_wTotal;
        
        
        % impact index
        CR_Avg.(expName).(aggTime).IwSnow  = C_wSnow./(C_wSnow+C_wSoil+C_wGW+C_wSurf);
        CR_Avg.(expName).(aggTime).IwSoil  = C_wSoil./(C_wSnow+C_wSoil+C_wGW+C_wSurf);
        CR_Avg.(expName).(aggTime).IwGW    = C_wGW./(C_wSnow+C_wSoil+C_wGW+C_wSurf);
        CR_Avg.(expName).(aggTime).IwSurf  = C_wSurf./(C_wSnow+C_wSoil+C_wGW+C_wSurf);
                
    end
 
end

[sname] = plotMaps_noCalc(expNames, varNames, CR, lat, lon, colLabel, spth_p)

colVars     = [rgb('RoyalBlue');rgb('SaddleBrown');rgb('DarkGreen');rgb('Purple')];
varNames    = {'IwSnow','IwSoil','IwGW','IwSurf'};
plotPieCR(expNames, varNames, CR_Avg, colVars, spth_p);


%% Time Series of TWS Components %%
% loop over experiments
for eN=1:numel(expNames)
    expName = expNames{eN};
    
    % get squeezed variables
    TS.fullD.(expName).wSnow = reshape(mod.(expName).wSnow, info.tem.helpers.sizes.nPix, info.tem.helpers.sizes.nTix);
    TS.fullD.(expName).wSoil = reshape(mod.(expName).wSoil(:,1,:), info.tem.helpers.sizes.nPix, info.tem.helpers.sizes.nTix) + reshape(mod.(expName).wSoil(:,2,:), info.tem.helpers.sizes.nPix, info.tem.helpers.sizes.nTix);
    TS.fullD.(expName).wDeep   = reshape(mod.(expName).wGW, info.tem.helpers.sizes.nPix, info.tem.helpers.sizes.nTix);
    TS.fullD.(expName).wSlow = reshape(mod.(expName).wSurf, info.tem.helpers.sizes.nPix, info.tem.helpers.sizes.nTix);
    TS.fullD.(expName).wTotal= reshape(mod.(expName).wTotal, info.tem.helpers.sizes.nPix, info.tem.helpers.sizes.nTix);
    
    TS.fullD.(expName).wOther        = TS.fullD.(expName).wDeep + TS.fullD.(expName).wSlow;
    TS.fullD.(expName).wTotal_noSnow = TS.fullD.(expName).wTotal - TS.fullD.(expName).wSnow;
    varNames = fieldnames(TS.fullD.(expName));
    
    % loop over the variables
    for vN=1:numel(varNames)
        varName = varNames{vN};
        % monthly values
        TS.full.(expName).(varName) = aggDay2Mon(TS.fullD.(expName).(varName),info.tem.model.time.sDate, info.tem.model.time.eDate);
        
        % only validObs?
        if withObs==1
            TS.full.(expName).(varName)(idxNaN) = NaN;
        end
 
        % detrend monthly data & anomalies
        TS.full.(expName).(varName)  = DetrendMatrix(TS.full.(expName).(varName));

        % calculate time anomalies of monthly data
        m = nanmean(TS.full.(expName).(varName),2);
        TS.ano.(expName).(varName) = TS.full.(expName).(varName) - repmat(m,1,length(xMonth));
        
        
        % calculate msc & iav for normal & for anomalies
        [tmp_msc,~, tmp_iav] = calcMSC(TS.ano.(expName).(varName),M);
        TS.msc_ano.(expName).(varName) = tmp_msc;
        TS.iav_ano.(expName).(varName) = tmp_iav;
        
        [tmp_msc,~, tmp_iav] = calcMSC(TS.full.(expName).(varName),M);
        TS.msc.(expName).(varName) = tmp_msc;
        TS.iav.(expName).(varName) = tmp_iav;
    end
end

% Plot the MSC with Impact Index
for eN=1:numel(expNames)
    expN = expNames{eN};
    dataCR.(expN).IwDeep     = CR.(expN).msc.IwGW;
    dataCR.(expN).IwSnow   = CR.(expN).msc.IwSnow;
    dataCR.(expN).IwSoil   = CR.(expN).msc.IwSoil;
    dataCR.(expN).IwSlow   = CR.(expN).msc.IwSurf;
end
CRnames     = {'IwSoil','IwDeep','IwSlow','IwSnow'}
varNames    = {'wTotal','wSoil','wDeep','wSlow','wSnow'};
colVars     = [rgb('Black');rgb('SaddleBrown');rgb('DarkGreen');rgb('Purple');rgb('RoyalBlue')];

[sname] = plotMSCvarsZones2(expNames, varNames, TS.msc_ano, [],dataCR,CRnames, zoneNames, KG_v, pix_a, [], colVars);
print(gcf,[spth_p 'ImpactIndex_' sname '.png'],'-dpng','-r300');


% Plot the IAV for zones 
for eN=1:numel(expNames)
    expN = expNames{eN};
    dataCR.(expN).IwDeep   = CR.(expN).iav.IwGW;
    dataCR.(expN).IwSnow   = CR.(expN).iav.IwSnow;
    dataCR.(expN).IwSoil   = CR.(expN).iav.IwSoil;
    dataCR.(expN).IwSlow   = CR.(expN).iav.IwSurf;
end

CRnames     = {'IwSoil','IwDeep','IwSlow','IwSnow'}
varNames    = {'wTotal','wSoil','wDeep','wSlow','wSnow'};
colVars     = [rgb('Black');rgb('SaddleBrown');rgb('DarkGreen');rgb('Purple');rgb('RoyalBlue')];
plotTSvarsZones2(expNames, varNames, TS.iav, [], dataCR, CRnames, zoneNames, KG_v, pix_a, xMonth, colVars, spth_p)

close all

%% -----------------------------------------------------------------
%% T over ET 
colXX  = [rgb('Blue');rgb('green')];
dataNames = expNames;

for eN = 1:numel(expNames)
    expN = expNames{eN};
    dataFull_ET.(expN).T_ET = mod.(expN).tranAct ./ mod.(expN).evapTotal;
    dataMSC_ET.(expN).T_ET  = calcMSC(dataFull_ET.(expN).T_ET, Md);
end

dataXX = {};
for op=1:numel(expNames)
    dataXX = [dataXX  dataMSC_ET.(expNames{op})];
end

[sname, sp] = plotMSCvarsZones_xExp_noObs(dataNames, dataXX, zoneNames, KG_v, pix_a, {'-', 'mm/d'}, colXX);
print(gcf,[spth  'MSC_ToverET.png'],'-dpng','-r300');

% Plot the mean T/ET Maps 
tmpLabel.colLim =  {[0.4 1]}; 
tmpLabel.colLimDiff  = [-0.2 0.2];

plotXMatrix_noObs(dataNames, dataXX, lat, lon,  tmpLabel, spth)

%% Plot the slow + fast + total runoff
clear dataXX

for eN = 1:numel(expNames)
    expN = expNames{eN};
    dataFull_Q.(expN).Qtotal = mod.(expN).roTotal;
    dataMSC_Q.(expN).Qtotal  = calcMSC(dataFull_Q.(expN).Qtotal, Md);
    dataFull_Q.(expN).Qslow  = mod.(expN).roSurfIndir;
    dataMSC_Q.(expN).Qslow   = calcMSC(dataFull_Q.(expN).Qslow, Md);
    dataFull_Q.(expN).Qfast  = mod.(expN).roSurfDir;
    dataMSC_Q.(expN).Qfast   = calcMSC(dataFull_Q.(expN).Qfast, Md);
end

dataXX = {};
for op=1:numel(expNames)
    dataXX = [dataXX  dataMSC_Q.(expNames{op})];
end

%plot the components of each experiment together
colXX  = [rgb('Black');rgb('Blue');rgb('green')];
vNames = fieldnames(dataMSC_Q.B);

obsMSC_Q.Q = dataObsMSC.Q;

plotMSCvarsZones2(expNames, vNames, dataMSC_Q, obsMSC_Q, [],[], zoneNames, KG_v, pix_a, {'mm/d', 'mm/d'}, colXX)
print(gcf,[spth  'MSC_allQ_Components.png'],'-dpng','-r300');

close all

%% --------------------------------------------------------------------- %%
%% PARAMTER VALUES
% load the parameters
for eN =1:numel(expNames)
    expName = expNames{eN};
    rDate   = rDates{eN};
    % load parameters
    params.(expName)  = load(['data/model_output/validCali_' expName '_validCali_' rDate '/modelOutput/paramValues.mat']);
end

% vegFrac
data2   = median(params.B.vegFrac,2);
expName = 'E_B_bL_RD4';
data    = median(params.E_B_bL_RD4.vegFrac,2);
figure, PlotMapGlobal_noInfo('Median Vegetation Fraction', ['VEG: median = ' num2str(round(mean(data(:),'omitnan'),2)) ' | B: ' num2str(round(mean(data2(:),'omitnan'),2)) ], [], lat, lon, data, [], othercolor('RdYlGn10'),[],1,1)
print(gcf,[spth 'VEG_vegFrac_median.png'],'-dpng','-r300');

% resulting smax2
data2   = params.B.p_wSoilBase_wAWC(:,2);
data    = params.E_B_bL_RD4.p_wSoilBase_wAWC(:,2);
figure, PlotMapGlobal_noInfo('Maximum Water Capacity - 2nd Soil Layer', ['VEG: median = ' num2str(round(median(data(:),'omitnan'),0)) 'mm | B: ' num2str(round(median(data2(:),'omitnan'),0)) 'mm'], 'mm', lat, lon, data, [], othercolor('BrBG10'),[],1,1)
print(gcf,[spth 'VEG_smax2_total.png'],'-dpng','-r300');


% scaled RD data
dataRD    = params.E_B_bL_RD4.p_wSoilBase_RD;
col       = othercolor('BrBG10');
colLim    = [0 prctile(dataRD(:),98)];
units     = 'mm';
nrows = 2;
ncols = 2;
% preps maps
land        = shaperead('landareas', 'UseGeoCoords', true);
rivers      = shaperead('worldrivers', 'UseGeoCoords', true);
geoRaRef    = georasterref('RasterSize', [180 360], 'RasterInterpretation', 'cells',  ...
    'LatitudeLimits', [-90 90], 'LongitudeLimits', [-180 180]);
Z1           = NaN(180, 360); % prepare mapgrid

figure; set(gcf, 'Position', [5 5 ncols*5 nrows*4]);
ha  = tight_subplot(nrows,ncols,[.02 .02],[.05 .05],[.05 .02]);
cnt = 1;
for rr=1:nrows
    for cc=1:ncols
        axes(ha(cnt))
        data        = dataRD(:,cnt);
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
        sm = surfm([-90 90], [-180 180], Z);
        geoshow(ha(cnt), rivers, 'Color', [.25 .25 .25]);
        t   =   title(['RD' num2str(cnt) ],'Fontsize', 7);
        colormap(ha(cnt),col);
        caxis(ha(cnt),colLim);
        if rr==nrows
            if cc==1
                cb  =   colorbar;
                cb1 =   ylabel(cb, units);
                set(cb,'FontSize',6, 'Position', [0.05 0.08 0.9 0.015],...
                    'Orientation', 'horizontal');
                set(cb1,'Position',[(colLim(2)+colLim(1))/2 -2 0]);
            end
        end
        cnt = cnt+1;
    end
end
%add the overall title name
an = annotation('textbox',[0 .5 1 .5],'String', ['VEG: Scaled Rooting Depth & Water Capacity Data'] ,'LineStyle','none', 'Fontsize', 10, 'Fontweight', 'b','fitboxtotext', 'on', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center');
print(gcf,[spth 'VEG_scaledRDdata.png'],'-dpng','-r300');

% Plot the effective alphaVeg
% a = p.tranDem.alphaVeg .* d.storedStates.vegFrac
% effAlphaVeg = tranDem / PET
timeX   = {info.tem.model.time.sDate,info.tem.model.time.eDate};
a       = infoIn.E_B_bL_RD4.info.tem.params.tranDem.alphaVeg .* squeeze(params.E_B_bL_RD4.vegFracMSC);

data = mean(a,2);
figure
PlotMapGlobal_noInfo('Mean Effective Alpha Coefficient (alphaVeg * vegFrac)','VEG',[],lat, lon,data,[],othercolor('RdYlGn10',250),[],1,1)
print(gcf,[spth  'VEG_effAlpha_mean.png'],'-dpng','-r300');

data = median(a,2);
PlotMapGlobal_noInfo('Median Effective Alpha Coefficient',['VEG | median = ' num2str(round(median(a(:)),2)) ],[],lat, lon,data,[],othercolor('RdYlGn10',250),[],1,1)
print(gcf,[spth expName 'VEG_effAlpha_median.png'],'-dpng','-r300');

% time series
plotTimeSeries2_gridstd(['Effective Alpha Coefficient (alphaVeg * vegFrac) | VEG' ],'d', timeX, {' ','effective alpha coefficient'},[], a)
print(gcf,[spth expName 'VEG_effAlpha_timeseries.png'],'-dpng','-r300');


