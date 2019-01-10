%Single Point Positioning Pragram which uses least sqaure method to estimate
%the GPS receiver's positionand clock offset at each epoch.This program's
%function includes:
%   1. Using PLSQ and pseudorange measurement to estimate the 3D single
%   point position for each epoch.Determine the position error by comparing
%   with given coordinates.Plot the "envelop" time series of the error
%   components.
%   2. Compute and plot the EDOP,NDOP,VDOP,HDOP,PDOP.
%   3. Compute and plot the residuals.
%
%Scripts written by Cheng Huang.

%Read observation and satellite files
[prn_s,t_s,x,y,z,xv,yv,zv] = readsat('Satellites.sat');
[n_obs,prn_o,t_o,pr,cp_l1,doppler_l1,cp_l2] = readobs('RemoteL1L2.obs');
XS = [x y z]; %Satellites Corrdinate (ECEF)

n_epoch = n_obs / 12; %total number of epoches
n_sat = zeros(n_epoch,1); %the number of satelltes on each epoch

for i = 1:n_epoch
    n_sat(i) = sum(prn_s((12*i-11):1:(12*i),:)~=0);
end

min_sat  = 4;

XR = zeros(n_epoch,3);
dT = zeros(n_epoch,1);
EDOP = zeros(n_epoch,1);
NDOP = zeros(n_epoch,1);
VDOP = zeros(n_epoch,1);
HDOP = zeros(n_epoch,1);
PDOP = zeros(n_epoch,1);
Cx_hat = zeros(n_epoch,3);
ELE = cell(n_epoch,1);
RES = cell(n_epoch,1);
PRN = cell(n_epoch,1);
n_sat_common = zeros(n_epoch,1);
pr_s = zeros(n_obs,1);
cp_s = zeros(n_obs,1);

%something wrong with sateliite file, prn_s, XS
for i = 1:n_epoch
    flag = find(prn_s(12*(i-1)+1:1:12*(i-1)+12));
    sat_flag = max(flag); %where satellite should end
    for j = 1:11
        if (prn_s(12*(i-1)+j) == 0 && prn_s(12*(i-1)+j+1) ~= 0)
            prn_s(12*(i-1)+j) = prn_s(12*(i-1)+j+1);
            prn_s(12*(i-1)+j+1) = 0;
            XS(12*(i-1)+j,:) = XS(12*(i-1)+j+1,:);
        end
    end
end

%Data Matching
for i = 1:n_epoch
    for k = 0:11
        for j = 0:11
            if (prn_s(12*i-11+k) == prn_o(12*i-11+j))
                pr_s(12*i-11+k) = pr(12*i-11+j);
                cp_s(12*i-11+k) = cp_l1(12*i-11+j);
            end
        end
    end
end

%===========================================================
%Iterative Least Square Formulation.
%---------------------------------------------------------------------------------------------------------
n_iteration = 0;
XR0 = [0 0 0]; %initial value for reveiver position
dT0 = 0;

for i = 1:n_epoch
    p1 = find(prn_o((12*(i-1)+1):1:12*i));
    for ii = 1:length(p1)
        p1(ii) = p1(ii) + 12*(i-1);
    end
    p2 = find(prn_s((12*(i-1)+1):1:12*i));
    for iii = 1:length(p2)
        p2(iii) = p2(iii) + 12*(i-1);
    end
    common_PRN = intersect(prn_o(p1),...
        prn_s(p2));
    PRN{i} = common_PRN;
    n_sat_common(i) = length(PRN{i});
    RES{i} = zeros(n_sat_common(i),1);
    index = zeros(n_sat_common(i),1);
    for j = 1:n_sat_common(i)
        for k = 1:12
            if prn_s(12*(i-1)+k) == PRN{i}(j)
                index(j) = 12*(i-1)+k;
            end
        end
    end
    if(n_sat(i) >= min_sat)
        while (n_iteration < 10)
            %update
            [XR(i,:),dT(i),EDOP(i),NDOP(i),VDOP(i),HDOP(i),PDOP(i),Cx_hat(i,:),RES{i}] = ...
                LS_adjustment(XR0, dT0,...
                XS(index,:), ...
                pr_s(index));
            
            delta_x = XR0 - XR(i,:);
            %exit judgement
            if (norm(delta_x(1:3)) < 1e-7)
                break;
            end
            
            XR0 = XR(i,:);
            dT0 = dT(i);
            n_iteration = n_iteration + 1;
        end
    end
    n_iteration = 0;
end
%===========================================================

%===========================================================
%Compare with the Reference Solution.
%---------------------------------------------------------------------------------------------------------
la_dms =  [51 15 31.11582];
lo_dms = [-114 06 01.76988];
[X_ref,Y_ref,Z_ref] = geodetic2ECEF(deg2rad(dms2degrees(la_dms)), ...
    deg2rad(dms2degrees(lo_dms)), ...
    1127.345, ...
    6378137, ...
    1/298.257223563); % WGS-84

x_ref = [X_ref;Y_ref;Z_ref];

diff_s = zeros(3,n_epoch);

for i = 1:n_epoch
    diff_s(:,i) = ECEF2ENU(x_ref,XR(i,:)');
end

%Determine the satellite elevation angle.
ele = zeros(n_obs,1);
count = 1;
for i = 1:n_obs
    if(i <= count*12 && pr_s(i) ~= 0)
        ele(i) = elevation(XR(count,:),XS(i,:));
    elseif(i > count*12 && pr_s(i) ~= 0)
        count = count + 1;
        ele(i) = elevation(XR(count,:),XS(i,:));
    end
end

for i = 1:n_epoch
    for j = 1:n_sat_common(i)
        for k = 1:12
            if prn_s(12*(i-1)+k) == PRN{i}(j)
                ELE{i}(j) = ele(12*(i-1)+k);
            end
        end
    end
end
%===========================================================

%creat a color map which contains 13 colors
cmap = linspecer(13);

tt = 1:1:n_epoch;
ttt = 1:1:n_obs;

%===========================================================
%Accuracy Comparison Plot.
%---------------------------------------------------------------------------------------------------------
figure(1)
h1 = subplot(3,1,1);
plot(tt,diff_s(1,:),'-r',tt,ones(n_epoch,1)*mean(diff_s(1,:)),'-b',tt,Cx_hat(:,1),'-k',tt,-Cx_hat(:,1),'-k','LineWidth',1);
set(gca,'fontsize',14)
grid on
xlabel('Epoch [s]','FontSize', 14)
ylabel('East Error [m]','FontSize', 14)
title(h1,'East Error','FontSize', 14)
legend(h1,{'True Accuracy','Mean of True Accuracy','Estimated Accuracy'})
h2 = subplot(3,1,2);
plot(tt,diff_s(2,:),'-r',tt,ones(n_epoch,1)*mean(diff_s(2,:)),'-b',tt,Cx_hat(:,2),'-k',tt,-Cx_hat(:,2),'-k','LineWidth',1);
grid on
xlabel('Epoch [s]','FontSize', 14)
ylabel('North Error [m]','FontSize', 14)
title(h2,'North Error','FontSize', 14)
h3 = subplot(3,1,3);
plot(tt,diff_s(3,:),'-r',tt,ones(n_epoch,1)*mean(diff_s(3,:)),'-b',tt,Cx_hat(:,3),'-k',tt,-Cx_hat(:,3),'-k','LineWidth',1 );
set(gca,'fontsize',14)
grid on
xlabel('Epoch [s]','FontSize', 14)
ylabel('Height Error [m]','FontSize', 14)
title(h3,'Height Error','FontSize', 14)
set(gcf,'units','points','position',[10,10,550,696])
set(gca,'fontsize',14)
hold on
%===========================================================

%===========================================================
%Satellite DOP and total number plot.
%---------------------------------------------------------------------------------------------------------
figure(2)
d1 = subplot(2,1,1);
plot(tt,EDOP,'b',tt,NDOP,'k',tt,VDOP,'r',tt,HDOP,'y',tt,PDOP,'g','LineWidth',2);
set(gca,'fontsize',14)
grid on
legend( 'EDOP','NDOP','VDOP','HDOP','PDOP');
xlabel('Epoch [s]','FontSize', 14)
title('DOP','FontSize', 14)
d2 = subplot(2,1,2);
plot(tt,n_sat,'.b','MarkerSize',10)
legend('Number of Satellites')
grid on
xlabel('Epoch [s]','FontSize', 14)
ylim([9 13])
title('Number of Satellites','FontSize', 14)
set(gcf,'units','points','position',[10,10,600,450])
set(gca,'fontsize',14)
hold on
%===========================================================

%===========================================================
%Residual Plot.
%---------------------------------------------------------------------------------------------------------
figure(3)
sat_type = intersect(prn_s,prn_o);
sat_type = sat_type(sat_type ~= 0); %exclude 0
n_sat_type = length(sat_type);

p7 = zeros(n_epoch,1);t7 = zeros(n_epoch,1);p8 = zeros(n_epoch,1);
t8 = zeros(n_epoch,1);p9 = zeros(n_epoch,1);t9 = zeros(n_epoch,1);
p11 = zeros(n_epoch,1);t11 = zeros(n_epoch,1);p15 = zeros(n_epoch,1);
t15 = zeros(n_epoch,1);p17 = zeros(n_epoch,1);t17 = zeros(n_epoch,1);
p18 = zeros(n_epoch,1);t18 = zeros(n_epoch,1);p19 = zeros(n_epoch,1);
t19 = zeros(n_epoch,1);p22 = zeros(n_epoch,1);t22 = zeros(n_epoch,1);
p24 = zeros(n_epoch,1);t24 = zeros(n_epoch,1);p26 = zeros(n_epoch,1);
t26 = zeros(n_epoch,1);p27 = zeros(n_epoch,1);t27 = zeros(n_epoch,1);
p28 = zeros(n_epoch,1);t28 = zeros(n_epoch,1);e7 = zeros(n_epoch,1);
e8 = zeros(n_epoch,1);e9 = zeros(n_epoch,1);e11 = zeros(n_epoch,1);
e15 = zeros(n_epoch,1);e17 = zeros(n_epoch,1);e18 = zeros(n_epoch,1);
e19 = zeros(n_epoch,1);e22 = zeros(n_epoch,1);e24 = zeros(n_epoch,1);
e26 = zeros(n_epoch,1);e27 = zeros(n_epoch,1);e28 = zeros(n_epoch,1);

for i = 1:n_epoch
    for j = 1:n_sat_common(i)
        if PRN{i}(j) == 7
            p7(i) = RES{i}(j);
            e7(i) = ELE{i}(j);
            t7(i) = i;
        elseif PRN{i}(j) == 8
            p8(i) = RES{i}(j);
            e8(i) = ELE{i}(j);
            t8(i) = i;
        elseif PRN{i}(j) == 9
            p9(i) = RES{i}(j);
            e9(i) = ELE{i}(j);
            t9(i) = i;
        elseif PRN{i}(j) == 11
            p11(i) = RES{i}(j);
            e11(i) = ELE{i}(j);
            t11(i) = i;
        elseif PRN{i}(j) == 15
            p15(i) = RES{i}(j);
            e15(i) = ELE{i}(j);
            t15(i) = i;
        elseif PRN{i}(j) == 17
            p17(i) = RES{i}(j);
            e17(i) = ELE{i}(j);
            t17(i) = i;
        elseif PRN{i}(j) == 18
            p18(i) = RES{i}(j);
            e18(i) = ELE{i}(j);
            t18(i) = i;
        elseif PRN{i}(j) == 19
            p19(i) = RES{i}(j);
            e19(i) = ELE{i}(j);
            t19(i) = i;
        elseif PRN{i}(j) == 22
            p22(i) = RES{i}(j);
            e22(i) = ELE{i}(j);
            t22(i) = i;
        elseif PRN{i}(j) == 24
            p24(i) = RES{i}(j);
            e24(i) = ELE{i}(j);
            t24(i) = i;
        elseif PRN{i}(j) == 26
            p26(i) = RES{i}(j);
            e26(i) = ELE{i}(j);
            t26(i) = i;
        elseif PRN{i}(j) == 27
            p27(i) = RES{i}(j);
            e27(i) = ELE{i}(j);
            t27(i) = i;
        elseif PRN{i}(j) == 28
            p28(i) = RES{i}(j);
            e28(i) = ELE{i}(j);
            t28(i) = i;
        end
    end
end
p7 = p7(p7 ~= 0);t7 = t7(t7 ~= 0);e7 = e7(e7 ~= 0);
p8 = p8(p8 ~= 0);t8 = t8(t8 ~= 0);e8 = e8(e8 ~= 0);
p9 = p9(p9 ~= 0);t9 = t9(t9 ~= 0);e9 = e9(e9 ~= 0);
p11 = p11(p11 ~= 0);t11 = t11(t11 ~= 0);e11 = e11(e11 ~= 0);
p15 = p15(p15 ~= 0);t15 = t15(t15 ~= 0);e15 = e15(e15 ~= 0);
% p17 = p17(p17 ~= 0);t17 = t17(t17 ~= 0);
p18 = p18(p18 ~= 0);t18 = t18(t18 ~= 0);e18 = e18(e18 ~= 0);
p19 = p19(p19 ~= 0);t19 = t19(t19 ~= 0);e19 = e19(e19 ~= 0);
p22 = p22(p22 ~= 0);t22 = t22(t22 ~= 0);e22 = e22(e22 ~= 0);
p24 = p24(p24 ~= 0);t24 = t24(t24 ~= 0);e24 = e24(e24 ~= 0);
p26 = p26(p26 ~= 0);t26 = t26(t26 ~= 0);e26 = e26(e26 ~= 0);
p27 = p27(p27 ~= 0);t27 = t27(t27 ~= 0);e27 = e27(e27 ~= 0);
t28 = t28(p28 ~= 0);p28 = p28(p28 ~= 0);e28 = e28(p28 ~= 0);

color_scheme = linspecer(13);
plot(t7,p7,'d','Marker','.','Color',color_scheme(1,:),'markers',12);hold on;
plot(t8,p8,'d','Marker','.','Color',color_scheme(2,:),'markers',12);hold on;
plot(t9,p9,'d','Marker','.','Color',color_scheme(3,:),'markers',12);hold on;
plot(t11,p11,'d','Marker','.','Color',color_scheme(4,:),'markers',12);hold on;
plot(t15,p15,'d','Marker','.','Color',color_scheme(5,:),'markers',12);hold on;
plot(t17,p17,'d','Marker','.','Color',color_scheme(6,:),'markers',12);hold on;
plot(t18,p18,'d','Marker','.','Color',color_scheme(7,:),'markers',12);hold on;
plot(t19,p19,'d','Marker','.','Color',color_scheme(8,:),'markers',12);hold on;
plot(t22,p22,'d','Marker','.','Color',color_scheme(9,:),'markers',12);hold on;
plot(t24,p24,'d','Marker','.','Color',color_scheme(10,:),'markers',12);hold on;
plot(t26,p26,'d','Marker','.','Color',color_scheme(11,:),'markers',12);hold on;
plot(t27,p27,'d','Marker','.','Color',color_scheme(12,:),'markers',12);hold on;
plot(t28,p28,'d','Marker','.','Color',color_scheme(13,:),'markers',12);hold on;

Legend = cell(13,1);
for iter=1:13
    Legend{iter}=strcat('PRN :G',num2str(sat_type(iter)));
end
legend(Legend)

grid on
xlabel('Epoch [s]','FontSize', 14)
ylabel('Residuals of Observations [m]','FontSize', 14)
title('Residuals of Observations','FontSize', 14)
set(gcf,'units','points','position',[10,10,660,360])
set(gca,'fontsize',14)
%===========================================================

%===========================================================
%Residual-Elevation Plot.
%---------------------------------------------------------------------------------------------------------
figure(4)
plot(e7,p7,'d','Marker','.','Color',color_scheme(1,:),'markers',12);hold on;
plot(e8,p8,'d','Marker','.','Color',color_scheme(2,:),'markers',12);hold on;
plot(e9,p9,'d','Marker','.','Color',color_scheme(3,:),'markers',12);hold on;
plot(e11,p11,'d','Marker','.','Color',color_scheme(4,:),'markers',12);hold on;
plot(e15,p15,'d','Marker','.','Color',color_scheme(5,:),'markers',12);hold on;
plot(e17,p17,'d','Marker','.','Color',color_scheme(6,:),'markers',12);hold on;
plot(e18,p18,'d','Marker','.','Color',color_scheme(7,:),'markers',12);hold on;
plot(e19,p19,'d','Marker','.','Color',color_scheme(8,:),'markers',12);hold on;
plot(e22,p22,'d','Marker','.','Color',color_scheme(9,:),'markers',12);hold on;
plot(e24,p24,'d','Marker','.','Color',color_scheme(10,:),'markers',12);hold on;
plot(e26,p26,'d','Marker','.','Color',color_scheme(11,:),'markers',12);hold on;
plot(e27,p27,'d','Marker','.','Color',color_scheme(12,:),'markers',12);hold on;
plot(e28,p28,'d','Marker','.','Color',color_scheme(13,:),'markers',12);hold on;
legend(Legend)

grid on
xlabel('Satellite Elevation Angle [deg]','FontSize', 14)
ylabel('Residuals of Observations [m]','FontSize', 14)
title('Residuals of Observations - Satellite Elevation Angle','FontSize', 14)
set(gcf,'units','points','position',[10,10,660,360])
set(gca,'fontsize',14)
%===========================================================
