function [cost] = calcMEF(obs, mod, unc)

tmp1    =   sum((obs-mod).^2./unc.^2, 'omitnan');
tmp2    =   sum((obs-mean(obs,'omitnan')).^2./unc.^2, 'omitnan');

cost    =   1-tmp1./tmp2;

end
