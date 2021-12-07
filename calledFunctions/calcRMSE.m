function [rmse] = calcRMSE(obs, mod)

tmp1    =   mean((obs-mod).^2, 'omitnan');

rmse    =   sqrt(tmp1);

end
