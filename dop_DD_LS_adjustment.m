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
%In LS solution, it only uses five satellites at every epoch which employs
%satellites with high DOP values.
%
%Scripts written by Cheng Huang.

function [XR,EDOP,NDOP,VDOP,HDOP,PDOP,Cx_hat,residual_obs,index] = dop_DD_LS_adjustment(XR0,X_ref,XS,pr_rem,pr_ref,pivot)

n = length(pr_rem); %number of obs
m = 3; %number of unknown parameters

if n<5
    disp('Number of Satellites is less than 5.');
    return;
end

%calculation of the distance vector
dis = zeros(n,1);
for i = 1:n
    dis(i) = sqrt(sum((XS(i,:)-XR0).^2,2)); %sum of each row, returns a column vector
end

dis_ref = zeros(n,1);
for i = 1:n
    dis_ref(i) = sqrt(sum((XS(i,:)-X_ref).^2,2)); %sum of each row, returns a column vector
end

XS_pivot = XS(pivot,:);
dis_pivot = dis(pivot);
dis_pivot_ref = dis_ref(pivot);
pr_rem_pivot = pr_rem(pivot);
pr_ref_pivot = pr_ref(pivot);

H1 = (XR0(1) - XS(:,1)) ./ dis - (XR0(1) - XS_pivot(1)) ./ dis_pivot; %column for X coordinates
H2 = (XR0(2) - XS(:,2)) ./ dis - (XR0(2) - XS_pivot(2)) ./ dis_pivot; %column for Y coordinate
H3 = (XR0(3) - XS(:,3)) ./ dis - (XR0(3) - XS_pivot(3)) ./ dis_pivot; %column for Z coordinate

for i = 1:n
    H1(i) = H1(i)*H1(i);
    H2(i) = H2(i)*H2(i);
    H3(i) = H3(i)*H3(i);
end

index = nchoosek((1:1:n),5); %all combinations of choosing 5 satellites from all the visible satellites in the epoch
n_combo = size(index,1);
sum_all = zeros(n_combo,1);
sum_x2 = zeros(n_combo,1);
sum_y2 = zeros(n_combo,1);
sum_z2 = zeros(n_combo,1);

for i = 1:n_combo
    sum_x2(i) = sum(H1(index(i,:)));
    sum_y2(i) = sum(H2(index(i,:)));
    sum_z2(i) = sum(H3(index(i,:)));
    
    sum_all(i) = 1 / sum_x2(i) + 1 / sum_y2(i) + 1 / sum_z2(i); %PDOP^2
end

[aa,k] = max(sum_all); %k is the index of the satellites that gives the highest PDOP values

XS = XS(index(k,:),:);
dis = dis(index(k,:));
dis_ref = dis_ref(index(k,:));
pr_rem = pr_rem(index(k,:));
pr_ref = pr_ref(index(k,:));

index = index(k,:);

%design matrix
H = [(XR0(1) - XS(:,1)) ./ dis - (XR0(1) - XS_pivot(1)) ./ dis_pivot, ... %column for X coordinate
	 (XR0(2) - XS(:,2)) ./ dis - (XR0(2) - XS_pivot(2)) ./ dis_pivot, ... %column for Y coordinate
	 (XR0(3) - XS(:,3)) ./ dis - (XR0(3) - XS_pivot(3)) ./ dis_pivot]; %column for Z coordinate

fx0 = (dis - dis_ref) - (dis_pivot - dis_pivot_ref); %fx0 = distance - cdt + iono-error + tropo_error

z = (pr_rem - pr_ref) - (pr_rem_pivot - pr_ref_pivot); %observation vector

% sigma_term = ones(n,1);
% R = diag(sigma_term);

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

%Computing DOP
Cx = sigma_hat*(H.'*H)^-1; 
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