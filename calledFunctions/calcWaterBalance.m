function [WB] = calcWaterBalance(inVar, stateVar, outVar)
%Checks Water Balance

% precision
preci = 1E-3;

% changes in storage
dS = diff(stateVar,1,2);

%subset I and O
I = inVar(:,2:end);
O = outVar(:,2:end);

WB = abs(I - O - dS);

WBCheck = WB > preci;

%info.checks.WBVioFrac = sum(sum(double(WBCheck)))/numel(WBCheck);
WBvioFrac = sum(sum(double(WBCheck)))/numel(WBCheck);

if WBvioFrac > 0
    %info.flags.WBalanceOK =0;
    mmsg=['Water balance not closed in ' num2str(WBvioFrac*100) ' % of cases'];
    warning(mmsg);
else
    disp('Water Balance is okay :)');
end


end % function
