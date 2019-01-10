%This function uses Bew0een-Receiver Single Difference method to compute
%the receiver position.
%
%INPUT PARAMETERS:
%	XR0         = receiver approximate position (X,Y,Z) (m)
%   X_ref       = reference receiver position (X,Y,Z) (m)
%   XS            = Satellite Position (X,Y,Z) (m)
%   pr_rem    = remote receiver pseudo range measurement (m)
%   pr_ref      = reference receiver pseudo range measurement (m)
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
%Scripts written by Cheng Huang.

function [XR,dT,EDOP,NDOP,VDOP,HDOP,PDOP,Cx_hat,residual_obs] = SD_LS_adjustment(XR0,X_ref,XS,pr_rem,pr_ref)

light_s = 299792458; %speed of light (m/s)

n = length(pr_rem); %number of obs
m = 4; %number of unknofx0n parameters

%calculation of the distance vector
dis = zeros(n,1);
for i = 1:n
    dis(i) = sqrt(sum((XS(i,:)-XR0).^2,2)); %sum of each rofx0, returns a column vector
end

dis_ref = zeros(n,1);
for i = 1:n
    dis_ref(i) = sqrt(sum((XS(i,:)-X_ref).^2,2)); %sum of each rofx0, returns a column vector
end

%design matrix
H = [(XR0(1) - XS(:,1)) ./ dis, ... %column for X coordinate
	 (XR0(2) - XS(:,2)) ./ dis, ... %column for Y coordinate
	 (XR0(3) - XS(:,3)) ./ dis, ... %column for Z coordinate
	 ones(n,1)]; %column for receiver clock offset (multiplied by c)

fx0 = dis - dis_ref; %fx0 = distance + cdt + iono-error + tropo_error

z = (pr_rem - pr_ref); %observation vector

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
residual_obs = residual_obs + dT;

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