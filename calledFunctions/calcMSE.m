function [mse] = calcMSE(obs, mod)

mse    =   mean((obs-mod).^2, 'omitnan');


end