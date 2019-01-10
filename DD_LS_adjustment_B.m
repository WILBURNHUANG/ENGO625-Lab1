function [XR,EDOP,NDOP,VDOP,HDOP,PDOP,Cx_hat,residual_obs] = DD_LS_adjustment_B(XR0,X_ref,XS,pr_rem,pr_ref,pivot)

n = length(pr_rem); %number of obs
m = 3; %number of unknown parameters

z = [pr_ref;pr_rem];

z1 = ones(n-1,1);
z2 = diag(ones(n-1,1));
B = [z1,-z2,-z1,z2];

zz = B*z;

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
H = [(XR0(1) - XS(2:n,1)) ./ dis(2:n) - (XR0(1) - XS(1,1)) ./ dis(1), ... %column for X coordinate
	 (XR0(2) - XS(2:n,2)) ./ dis(2:n) - (XR0(2) - XS(1,2)) ./ dis(1), ... %column for Y coordinate
	 (XR0(3) - XS(2:n,3)) ./ dis(2:n) - (XR0(3) - XS(1,3)) ./ dis(1)]; %column for Z coordinate

fx0 = (dis(2:n) - dis_ref(2:n)) - (dis(1) - dis_ref(1)); %fx0 = distance - cdt + iono-error + tropo_error

% z = (pr_rem - pr_ref) - (pr_rem(pivot) - pr_ref(pivot)); %observation vector

sigma_term = ones(n-1,1);
P = diag(sigma_term);

N = H'*P*H; %normal matrix

%LS formulation
x = N \ (H'*P*(zz - fx0)); 

z_hat = H*x + fx0;
r_hat = zz - z_hat;

Cx_hat = zeros(1,3);

%output parameters calculation
sigma_hat = (r_hat'*r_hat) / (n-1-m);
residual_obs = [0;zz - fx0 - H*x];
XR = XR0 + x(1:3).';

%Computing DOP
Cx = sigma_hat*(H'*(B*eye(2*n)*B')^-1*H)^-1; 
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