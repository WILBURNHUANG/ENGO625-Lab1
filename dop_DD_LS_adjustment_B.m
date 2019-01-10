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

function [XR,EDOP,NDOP,VDOP,HDOP,PDOP,Cx_hat,residual_obs,index] = dop_DD_LS_adjustment_B(XR0,X_ref,XS,pr_rem,pr_ref,pivot)

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

H1 = (XR0(1) - XS(2:n,1)) ./ dis(2:n) - (XR0(1) - XS(1,1)) ./ dis(1); %column for X coordinates
H2 = (XR0(2) - XS(2:n,2)) ./ dis(2:n) - (XR0(2) - XS(1,2)) ./ dis(1); %column for Y coordinate
H3 = (XR0(3) - XS(2:n,3)) ./ dis(2:n) - (XR0(3) - XS(1,3)) ./ dis(1); %column for Z coordinate

for i = 1:n-1
    H1(i) = H1(i)*H1(i);
    H2(i) = H2(i)*H2(i);
    H3(i) = H3(i)*H3(i);
end

index = nchoosek((1:1:n-1),4); %all combinations of choosing 5 satellites from all the visible satellites in the epoch
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

XS = [XS(1,:);XS(index(k,:),:)];
dis =  [dis(1);dis(index(k,:))];
dis_ref = [dis_ref(1);dis_ref(index(k,:))];
pr_rem = [pr_rem(1);pr_rem(index(k,:))];
pr_ref = [pr_ref(1);pr_ref(index(k,:))];

index = [1,index(k,:)];

z = [pr_ref;pr_rem];

z1 = ones(4,1);
z2 = diag(ones(4,1));
B = [z1,-z2,-z1,z2];

zz = B*z;

%design matrix
H = [(XR0(1) - XS(2:5,1)) ./ dis(2:5) - (XR0(1) - XS(1,1)) ./ dis(1), ... %column for X coordinate
	 (XR0(2) - XS(2:5,2)) ./ dis(2:5) - (XR0(2) - XS(1,2)) ./ dis(1), ... %column for Y coordinate
	 (XR0(3) - XS(2:5,3)) ./ dis(2:5) - (XR0(3) - XS(1,3)) ./ dis(1)]; %column for Z coordinate

fx0 = (dis(2:5) - dis_ref(2:5)) - (dis(1) - dis_ref(1)); %fx0 = distance - cdt + iono-error + tropo_error

% sigma_term = ones(n,1);
% R = diag(sigma_term);

N = H'*H; %normal matrix

%LS formulation
x = N \ (H'*(zz - fx0)); 

z_hat = H*x + fx0;
r_hat = zz - z_hat;

Cx_hat = zeros(1,3);

%output parameters calculation
sigma_hat = (r_hat'*r_hat) / (n-m);
residual_obs = [0;zz - fx0 - H*x];
XR = XR0 + x(1:3).';

%Computing DOP
Cx = sigma_hat*(H'*(B*eye(10)*B')^-1*H)^-1; 
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