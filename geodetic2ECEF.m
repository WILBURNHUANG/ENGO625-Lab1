function [X,Y,Z] = geodetic2ECEF (phi, lam, h, a, f)
e = sqrt(1-(1-f)^2);
e2 = 1 - (1 - f)^2;

N = a ./ sqrt(1 - e2 * sin(phi).^2);

X = (N + h) .* cos(lam) .* cos(phi);
Y = (N + h) .* sin(lam) .* cos(phi);
Z = (N * (1 - e2) + h) .* sin(phi);
