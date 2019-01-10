function [El] = elevation(XR, XS)
[phi, lam] = ECEF2geodetic(XR(1), XR(2), XR(3));

%rotation matrix from global to local reference system
R = [-sin(lam) cos(lam) 0;
    -sin(phi)*cos(lam) -sin(phi)*sin(lam) cos(phi);
    +cos(phi)*cos(lam) +cos(phi)*sin(lam) sin(phi)];

local = R*(XS-XR)';
E = local(1);
N = local(2);
U = local(3);
dis = sqrt(E^2+N^2);

if dis < 1.e-20
   El = 90;
else
   El = rad2deg(atan2(U,dis));
end

