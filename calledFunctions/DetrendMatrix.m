
function [MyData,ltr]=DetrendMatrix(MyData)
sz=size(MyData);
ltr=NaN(sz(1),1); % linear trend for each pixel (in rows)
yrs=1:sz(2); % number of timesteps (in cols)

%col
for row=1:sz(1) % iterate over pixel
%     p = polyfit(yrs,MyData(row,:),1); % polynomial fit of degree 1 (= linear) of input data
%     tr=yrs.*p(1)+p(2); % calculate linear trend from polyfit coefficients for each time step
%     MyData(row,:)=MyData(row,:)-tr; % remove trend from input data
%     ltr(row,1)=p(1); % save trend coefficient
%     %p(1)
    
        % without NaNs
    data_tmp  = MyData(row,:);
    Good      = isnan(yrs) + isnan(data_tmp);
    p         = polyfit(yrs(Good==0),data_tmp(Good==0),1);
    tr        = yrs(Good==0).*p(1)+p(2); % calculate linear trend from polyfit coefficients for each time step
    data_tmp(Good==0)  = data_tmp(Good==0)-tr; % remove trend from input data
%%%% what to do with NaN points????? %%%%%%%
    data_tmp(Good ~=0) = NaN;
    MyData(row,:) = data_tmp;
    ltr(row,1)    = p(1); % save trend coefficient

end



end

% % original
% sz=size(MyData);
% ltr=NaN(1,sz(2));
% yrs=1:sz(1);
% 
% for row=1:sz(2)
%     
%     p = polyfit(yrs',squeeze(MyData(:,row)),1);
%     tr=yrs.*p(1)+p(2);
%     MyData(:,row)=squeeze(MyData(:,row))-tr';
%     ltr(row)=p(1);
% end
