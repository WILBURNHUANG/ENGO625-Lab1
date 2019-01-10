%This function uses Parametric Least Square (PLSQ) method and pesudorange
%measurenment to estimate the 3D single point position for each epoch.
%In LS solution, it only uses five satellites at every epoch ehich employs
%satellites eith high DOP values.
%
%Scripts written by Cheng Huang.

function [XR,dT,EDOP,NDOP,VDOP,HDOP,PDOP,Cx_hat,residual_obs,index] = dop_LS(XR0,dT0,XS,pr)

light_s = 299792458; %speed of light (m/s)

n = length(pr); %number of obs
m = 4; %number of unknofx0n parameters

if n<5
    disp('Number of Satellites is less than 5.');
    return;
end

%calculation of the distance vector
distance = zeros(n,1);
for i = 1:n
    distance(i) = sqrt(sum((XS(i,:)-XR0).^2,2)); %sum of each rofx0, returns a column vector
end

H1 = (XR0(1) - XS(:,1)) ./ distance; %column for X coordinates
H2 = (XR0(2) - XS(:,2)) ./ distance; %column for Y coordinate
H3 = (XR0(3) - XS(:,3)) ./ distance; %column for Z coordinate

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
distance = distance(index(k,:));
pr = pr(index(k,:));

index = index(k,:);

%design matrix
H = [(XR0(1) - XS(:,1)) ./ distance, ... %column for X coordinate
	 (XR0(2) - XS(:,2)) ./ distance, ... %column for Y coordinate
	 (XR0(3) - XS(:,3)) ./ distance, ... %column for Z coordinate
	 ones(5,1)]; %column for receiver clock offset (multiplied by c)

fx0 = distance - light_s*dT0; %fx0 = distance - cdt + iono-error + tropo_error

z = pr; %observation vector

N = H'*H; %normal matrix

%Computing DOP
Cx = (H.'*H)^-1; %Cofactor Matrix
Cx = Cx(1:3,1:3); %the fourth term is about time precision
Q_enu = ECEF2ENU_Cov(Cx, XR0(1,:)');
EDOP = sqrt(Q_enu(1,1));
NDOP = sqrt(Q_enu(2,2));
VDOP = sqrt(Q_enu(3,3));
HDOP = sqrt(Q_enu(1,1) + Q_enu(2,2));
PDOP = sqrt(Q_enu(1,1) + Q_enu(2,2) + Q_enu(3,3));

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

Cx_hat(1) = sqrt(sigma_hat*Q_enu(1,1));
Cx_hat(2) = sqrt(sigma_hat*Q_enu(2,2));
Cx_hat(3) = sqrt(sigma_hat*Q_enu(3,3));
end