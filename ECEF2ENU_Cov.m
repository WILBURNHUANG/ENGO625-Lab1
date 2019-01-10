function [Cyy] = ECEF2ENU_Cov(Cxx, X)
%initialize new covariance matrices
Cyy = zeros(size(Cxx));

for i = 1 : size(X,2)

    %geodetic coordinates
    [phi, lam] = ECEF2geodetic(X(1,i), X(2,i), X(3,i));

    %rotation matrix from global to local reference system
    R = [-sin(lam) cos(lam) 0;
         -sin(phi)*cos(lam) -sin(phi)*sin(lam) cos(phi);
         +cos(phi)*cos(lam) +cos(phi)*sin(lam) sin(phi)];

    %covariance propagation
    Cyy(:,:,i) = R * Cxx(:,:,i) * R';
end