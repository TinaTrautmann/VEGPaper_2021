function [KGE, r_p, alpha, beta] = calcKGE(obs, mod)

% calculates the Kling Gupta Efficency (Gupta et al. 2009)
% output:
%       KGE     - Kling Gupta Efficency
%       r_p     - Pearson correlation
%       alpha   - variability error
%       beta    - bias error
% input:
%       obs     - observations (rows = npix, cols = ntix)
%       mod     - simulations (rows = npix, cols = ntix)


%% use consistent data points
v_obs   = find(~isnan(obs));

if isempty(v_obs)
    warning('no v_obs!');
    r_p = NaN;
    alpha = NaN;
    beta = NaN;
    KGE = NaN;
else
    
    
    
    %% calculations
    r_p     =   corr(obs(v_obs), mod(v_obs), 'rows', 'complete');
    
    alpha   =   std(mod(v_obs)) ./ std(obs(v_obs));
    
    beta    =   mean(mod(v_obs)) ./ mean(obs(v_obs));
    
    % Euclidian distance
    ED      =   sqrt((r_p - 1).^2 + (alpha  - 1).^2 + (beta  - 1).^2);
    
    % Kling Gupta Efficency
    KGE     =   1 - ED;
    
end

end


