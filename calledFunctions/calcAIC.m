function [AIC, BIC] = calcAIC(obs, mod, n_params)

% obs       - observations, all in one column
% mod       - different columns = different models
% n_params  - number of params of each model, also per column

% AIC(1,:)  - AIC index
% AIC(2,:)  - delta AIC (difference to minimum)

% number of data points
n_obs       = find(~isnan(obs));
n_tmp       = length(n_obs);

n_exp = size(mod,2);

AIC = NaN(2,n_exp);
BIC = NaN(2,n_exp);

for nE = 1:n_exp
     p_tmp       = n_params(1,nE);
     mod_tmp     = mod(:,nE);
     MSE_tmp     = calcMSE(obs(:), mod_tmp);
     
    [AIC_tmp,BIC_tmp]   = MyAICandBIC(MSE_tmp,n_exp,p_tmp);
    
    AIC(1,nE) = AIC_tmp;
    BIC(1,nE) = BIC_tmp;
end

% deltaAIC (difference to min)
aic_min = min(AIC(1,:),[],2);
bic_min = min(BIC(1,:),[],2);

for nE = 1:n_exp
    AIC(2,nE) = AIC(1,nE) - aic_min;
    BIC(2,nE) = BIC(1,nE) - bic_min;
end

end