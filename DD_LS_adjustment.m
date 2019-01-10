%This function uses Between-Receiver Double Difference method to compute
%the receiver position.
%
%INPUT PARAMETERS:
%	XR0         = receiver approximate position (X,Y,Z) (m)
%   X_ref       = reference receiver position (X,Y,Z) (m)
%   XS            = Satellite Position (X,Y,Z) (m)
%   pr_rem    = remote receiver pseudo range measurement (m)
%   pr_ref      = reference receiver pseudo range measurement (m)
%   pivot       = index for the satellite with highest elevation angle
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

function [XR,EDOP,NDOP,VDOP,HDOP,PDOP,Cx_hat,residual_obs] = DD_LS_adjustment(XR0,X_ref,XS,pr_rem,pr_ref,pivot)

n = length(pr_rem); %number of obs
m = 3; %number of unknown parameters

%calculation of the distance vector
dis = zeros(n,1);
for i = 1:n
    dis(i) = sqrt(sum((XS(i,:)-XR0).^2,2)); %sum of each row, returns a column vector
end

dis_ref = zeros(n,1);
for i = 1:n
    dis_ref(i) = sqrt(sum((XS(i,:)-X_ref).^2,2)); %sum of each row, returns a column vector
end

%design matrix
H = [(XR0(1) - XS(:,1)) ./ dis - (XR0(1) - XS(pivot,1)) ./ dis(pivot), ... %column for X coordinate
	 (XR0(2) - XS(:,2)) ./ dis - (XR0(2) - XS(pivot,2)) ./ dis(pivot), ... %column for Y coordinate
	 (XR0(3) - XS(:,3)) ./ dis - (XR0(3) - XS(pivot,3)) ./ dis(pivot)]; %column for Z coordinate

fx0 = (dis - dis_ref) - (dis(pivot) - dis_ref(pivot)); %fx0 = distance - cdt + iono-error + tropo_error

z = (pr_rem - pr_ref) - (pr_rem(pivot) - pr_ref(pivot)); %observation vector

sigma_term = ones(n,1);
R = diag(sigma_term);

N = H'*R*H; %normal matrix

%LS formulation
x = N \ (H'*R*(z - fx0)); 

z_hat = H*x + fx0;
r_hat = z - z_hat;

Cx_hat = zeros(1,3);

%output parameters calculation
sigma_hat = (r_hat'*r_hat) / (n-m);
residual_obs = z - fx0 - H*x;
XR = XR0 + x(1:3).';

%Computing DOP
Cx = sigma_hat*(H.'*R*H)^-1; 
Cx = Cx(1:3,1:3); %the fourth term is about time precision
Cx_enu = ECEF2ENU_Cov(Cx, XR(1,:)');

Q_x = (H.'*H); %Cofactor Matrix
Q_enu = ECEF2ENU_Cov(Q_x, XR(1,:)');
EDOP = sqrt(Q_enu(1,1));
NDOP = sqrt(Q_enu(2,2));
VDOP = sqrt(Q_enu(3,3));
HDOP = sqrt(Q_enu(1,1) + Q_enu(2,2));
PDOP = sqrt(Q_enu(1,1) + Q_enu(2,2) + Q_enu(3,3));

Cx_hat(1) = sqrt(Cx_enu(1,1));
Cx_hat(2) = sqrt(Cx_enu(2,2));
Cx_hat(3) = sqrt(Cx_enu(3,3));
end