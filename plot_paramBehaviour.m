%% set save pth
spth = [pwd '/plots/parameter_behaviour/'];
if ~exist(spth, 'dir'), mkdir(spth),end

%% --- p.berg (fixed. p.smaxVeg)
p.smaxVeg       = 200;
p.berg          = [0:0.1:10]
    
WBP             = repmat([0:0.1:50]',1,length(p.berg));
SoilTotal       = NaN(size(WBP));
SoilTotal(1,:)  = 0;
InfExFrac       = NaN(size(WBP));
InSoil          = NaN(size(WBP));
QoverFlow       = NaN(size(WBP));

for bb=1:length(p.berg)
    for ii=1:size(WBP,1)
        InfExFrac(ii,bb)   = min(exp(p.berg(bb) .* log(SoilTotal(ii,bb)  ./ p.smaxVeg)),1);
        QoverFlow(ii,bb)   = WBP(ii,bb) .* InfExFrac(ii,bb);
        
        InSoil(ii,bb)        = WBP(ii,bb) .* (1-InfExFrac(ii,bb));
        SoilTotal(ii+1,bb)     = SoilTotal(ii,bb) + InSoil(ii,bb);
    end
end

Q_WBP  = QoverFlow ./ WBP;
SM_max = SoilTotal ./ p.smaxVeg;
SM_max = SM_max(1:end-1,:);

col = jet(length(p.berg));

nameT = ['fix smax at ' num2str(p.smaxVeg) ' mm'];
dataX = SM_max;
dataY = Q_WBP;
nameX = 'wSoil/smax';
nameY = 'Qoverflow/WBP';

figure
for bb=1:length(p.berg)
    plot(dataX(:,bb), dataY(:,bb), 'color', col(bb,:))
    hold on
end
cb = colorbar, colormap(jet(length(p.berg))), caxis([min(p.berg) max(p.berg)]), cb.Label.String = 'p.berg';
ylabel(nameY), xlabel(nameX),title(nameT), set(gca, 'fontsize', 9)
print(gcf, [spth nameT '.png'], '-dpng', '-r300');

%% --- p.smaxVeg (fixed p.berg)
p.smaxVeg       = [0:0.1:500];
p.berg          = 2;
    
WBP             = repmat([0:0.1:50]',1,length(p.smaxVeg));
SoilTotal       = NaN(size(WBP));
SoilTotal(1,:)  = 0;
InfExFrac       = NaN(size(WBP));
InSoil          = NaN(size(WBP));
QoverFlow       = NaN(size(WBP));

for bb=1:length(p.smaxVeg)
    for ii=1:size(WBP,1)
        InfExFrac(ii,bb)   = min(exp(p.berg .* log(SoilTotal(ii,bb)  ./ p.smaxVeg(bb))),1);
        QoverFlow(ii,bb)   = WBP(ii,bb) .* InfExFrac(ii,bb);
        
        InSoil(ii,bb)        = WBP(ii,bb) .* (1-InfExFrac(ii,bb));
        SoilTotal(ii+1,bb)     = SoilTotal(ii,bb) + InSoil(ii,bb);
    end
end

Q_WBP  = QoverFlow ./ WBP;
SM_max = SoilTotal ./ p.smaxVeg;
SM_max = SM_max(1:end-1,:);

col = jet(length(p.smaxVeg));

nameT = ['fix berg at ' num2str(p.berg) ' mm'];
dataX = SoilTotal(1:end-1,:);
dataY = Q_WBP;
nameX = 'wSoil';
nameY = 'Qoverflow/WBP';
nameC = 'p.smaxVeg'

figure
for bb=1:length(p.smaxVeg)
    plot(dataX(:,bb), dataY(:,bb), 'color', col(bb,:))
    hold on
end
cb = colorbar, colormap(jet(length(p.smaxVeg))), caxis([min(p.smaxVeg) max(p.smaxVeg)]), cb.Label.String = nameC;
ylabel(nameY), xlabel(nameX),title(nameT), set(gca, 'fontsize', 9)
print(gcf, [spth nameT '.png'], '-dpng', '-r300');

