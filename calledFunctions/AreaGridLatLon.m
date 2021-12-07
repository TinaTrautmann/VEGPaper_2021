% Created by Naixin Fan. The algrithm is from Martin Jung. 
% Thank Uli for sharing his script.
% This is script is for calculating the area of grid cells. 

function area_grid = AreaGridLatLon(lats,lons,res)


res_lat     = res(1);
res_lon     = res(2);

% ER          =   6371.004; %Earth radius (km)
ER          =   6378160; %Earth radius (m)

cols        =   size(lons,1);

londel=abs(res_lon);
lats1 = lats - res_lat/2;
lats2 = lats + res_lat/2;
areavec=(pi/180)*ER^2*abs(sin(lats1.*pi/180)-sin(lats2.*pi/180))*londel;

area_grid       = areavec * ones(1, cols);
end
