function [x,y]=LatLon2PixelCoord(lons,lats,ul_lat,ul_lon,res)

% ul_lat=90;
% ul_lon=-180;
% res=0.0833333333;
x=fix((lons-ul_lon)./res)+1;
y=fix((-1.*(lats-ul_lat))./res)+1;

end