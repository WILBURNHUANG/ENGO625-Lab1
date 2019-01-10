%This function uses Parametric Least Square (PLSQ) method and pesudorange
%measurenment to estimate the 3D single point position for each epoch.
%
%INPUT PARAMETERS:
%	XR0	= receiver approximate position (X,Y,Z) (m)
%   XS      = Satellite Position (X,Y,Z) (m)
%   pr      = pseudo range measurement (m)
%
%OUTPUT :  
%   XR      = receiver position after estimation (X,Y,Z) (m)
%   dT      = receiver clock offset
%   EDOP
%   NDOP
%   VDOP
%   HDOP
%   PDOP
%   sigma_hat   = posterior sigma
%   residual_obs    = residuals of all input observation
%
%Scripts fx0ritten by Cheng Huang.

function [XR,dT,EDOP,NDOP,VDOP,HDOP,PDOP,Cx_hat,residual_obs] = LS_adjustment(XR0,dT0,XS,pr)

light_s = 299792458; %speed of light (m/s)

n = length(pr); %number of obs
m = 4; %number of unknofx0n parameters

%calculation of the distance vector
distance = zeros(n,1);
for i = 1:n
    distance(i) = sqrt(sum((XS(i,:)-XR0).^2,2)); %sum of each rofx0, returns a column vector
end

%design matrix
H = [(XR0(1) - XS(:,1)) ./ distance, ... %column for X coordinate
	 (XR0(2) - XS(:,2)) ./ distance, ... %column for Y coordinate
	 (XR0(3) - XS(:,3)) ./ distance, ... %column for Z coordinate
	 ones(n,1)]; %column for receiver clock offset (multiplied by c)

fx0 = distance - light_s*dT0; %fx0 = distance - cdt + iono-error + tropo_error

z = pr; %observation vector

N = H'*H; %normal matrix

%LS formulation
x = N \ (H'*(z - fx0)); 

z_hat = H*x + fx0;
r_hat = z - z_hat;

Cx_hat = zeros(1,3);

%output parameters calculation
sigma_hat = (r_hat'*r_hat) / (n-m);
residual_obs = z - fx0 - H*x;
XR = XR0 + x(1:3).';
dT = x(4) / light_s;

%Computing DOP
Cx = (H.'*H)^-1; %Cofactor Matrix
Cx = Cx(1:3,1:3); %the fourth term is about time precision
Q_enu = ECEF2ENU_Cov(Cx, XR(1,:)');
EDOP = sqrt(Q_enu(1,1));
NDOP = sqrt(Q_enu(2,2));
VDOP = sqrt(Q_enu(3,3));
HDOP = sqrt(Q_enu(1,1) + Q_enu(2,2));
PDOP = sqrt(Q_enu(1,1) + Q_enu(2,2) + Q_enu(3,3));

Cx_hat(1) = sqrt(sigma_hat*Q_enu(1,1));
Cx_hat(2) = sqrt(sigma_hat*Q_enu(2,2));
Cx_hat(3) = sqrt(sigma_hat*Q_enu(3,3));
end