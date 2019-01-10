%Between-receiver-single-differencing mode to compute the position of
%remote receiver. There are a remote receiver and a reference receiver that
%keep track of the same satellite at the same epoch. By subtracting the
%pesudorange at reference station from that at remote, orbital errors and
%atmospheric errors are reduced, satellite clock error are eliminated.
%
%Scripts written by Cheng Huang.

%Read observation and satellite files
[prn_s,t_s,x,y,z,xv,yv,zv] = readsat('Satellites.sat');
[n_obs,prn_rem,t_rem,pr_rem,cp_l1_rem,doppler_l1_rem,cp_l2_rem] = readobs('RemoteL1L2.obs');
[n_obs_ref,prn_ref,t_ref,pr_ref,cp_l1_ref,doppler_l1_ref,cp_l2_ref] = readobs('BaseL1L2.obs');

XS = [x y z]; %Satellites Corrdinate (ECEF)
n_epoch = n_obs / 12;
n_epoch_ref = n_obs_ref / 12;
n_sat = zeros(n_epoch,1);
n_sat_ref = zeros(n_epoch_ref,1);

for i = 1:n_epoch
    n_sat(i) = sum(prn_s((12*i-11):1:(12*i),:)~=0);
end

for i = 1:n_epoch_ref
    n_sat_ref(i) = sum(prn_s((12*i-11):1:(12*i),:)~=0);
end

min_sat  = 4;

%residual_obs = zeros(n_obs,1);
pr_s = zeros(n_obs,1);
pr_s_ref = zeros(n_obs_ref,1);

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
            if (prn_s(12*i-11+k) == prn_rem(12*i-11+j))
                pr_s(12*i-11+k) = pr_rem(12*i-11+j);
            end
            if (prn_s(12*i-11+k) == prn_ref(12*i-11+j))
                pr_s_ref(12*i-11+k) = pr_ref(12*i-11+j);
            end
        end
    end
end

%===========================================================
%SPP for  Reference Receiver.
%---------------------------------------------------------------------------------------------------------
X_ref = zeros(n_epoch_ref,3);
n_iteration = 0;
XR0_ref = [0 0 0]; %initial value for reveiver position
dT0_ref = 0;

for i = 1:n_epoch_ref
    index_ref = find(pr_s_ref((12*i-11):1:(12*i - 12 + n_sat_ref(i))));
    for t = 1:length(index_ref)
        index_ref(t) = index_ref(t) + 12*(i-1);
    end
    if(n_sat_ref(i) >= min_sat)
        while (n_iteration < 10)
            %update
            [X_ref(i,:), ~ , ~  , ~ , ~ , ~ ] = LS_adjustment(XR0_ref, dT0_ref,XS(index_ref,:), pr_s_ref(index_ref));
            delta_x = XR0_ref - X_ref(i,:);
            %exit judgement
            if (norm(delta_x(1:3)) < 1e-7)
                break;
            end
            
            XR0_ref = X_ref(i,:);
            n_iteration = n_iteration + 1;
        end
    end
    n_iteration = 0;
end
%===========================================================

%===========================================================
%Synchronizing remote receiver and reference receiver.
%---------------------------------------------------------------------------------------------------------
index_start_t = 0;
index_end_t = 0;
if t_ref(1) <= t_rem(1) && t_ref(n_obs_ref-11) >= t_rem(n_obs-11)
    %reference receiver time range can cover the remote receiver
    for i = 1:n_obs_ref
        if (t_ref(i) == t_rem(1))
            index_start_t = i;
            break;
        end
    end
    
    for k = n_obs_ref:-1:1
        if (t_ref(k) == t_rem(n_obs-11))
            index_end_t = k;
            if index_end_t / 12 ~= 0
                index_end_t = (floor(index_end_t/12) + 1)*12;
            end
            break;
        end
    end
    
    pr_s_ref = pr_s_ref(index_start_t:1:index_end_t);
    start_epoch = (index_start_t-1)/12 + 1;
    end_epoch = index_end_t/12;
    X_ref = X_ref(start_epoch:1:end_epoch,:);
    XS = XS(index_start_t:1:index_end_t,:);
    prn_ref = prn_ref(index_start_t:1:index_end_t);
end

if (t_ref(1) > t_rem(1) && t_ref(n_obs_ref-11) < t_rem(n_obs-11))
    %remote receiver time range can cover the reference receiver
    for i = 1:n_obs
        if (t_rem(i) == t_ref(1))
            index_start_t = i;
            break;
        end
    end
    
    for k = n_obs:-1:1
        if (t_rem(k) == t_ref(n_obs_ref))
            index_end_t = k;
            break;
        end
    end
    start_epoch = (index_start_t-1)/12 + 1;
    end_epoch = index_end_t/12;
    pr_s = pr_s(index_start_t:1:index_end_t);
    n_epoch = n_epoch_ref;
    n_obs = n_obs_ref;
    XS = XS(index_start_t:1:index_end_t,:);
    prn_rem = prn_ref(index_start_t:1:index_end_t);
end

if (t_ref(1) <= t_rem(1) && t_ref(n_obs_ref-11) < t_rem(n_obs-11))
    %cover the start not the end
    for i = 1:n_obs_ref
        if (t_ref(i) == t_rem(1))
            index_start_t = i;
            break;
        end
    end
    
    index_end_t = n_obs_ref;
    
    pr_s_ref = pr_s_ref(index_start_t:1:index_end_t);
    start_epoch = (index_start_t-1)/12 + 1;
    end_epoch = index_end_t/12;
    X_ref = X_ref(start_epoch:1:end_epoch,:);
    XS = XS(index_start_t:1:index_end_t,:);
    prn_ref = prn_ref(index_start_t:1:index_end_t);
    prn_rem = prn_rem(index_start_t:1:index_end_t);
end

if (t_ref(1) > t_rem(1) && t_ref(n_obs_ref-11) >= t_rem(n_obs-11))
    %cover the end not the start
    index_start_t = 1;
    
    for k = n_obs_ref:-1:1
        if (t_ref(k) == t_rem(n_obs))
            index_end_t = k;
            break;
        end
    end
    
    pr_s_ref = pr_s_ref(index_start_t:1:index_end_t);
    start_epoch = (index_start_t-1)/12 + 1;
    end_epoch = index_end_t/12;
    X_ref = X_ref(start_epoch:1:end_epoch,:);
    XS = XS(index_start_t:1:index_end_t,:);
    prn_ref = prn_ref(index_start_t:1:index_end_t);
    prn_rem = prn_rem(index_start_t:1:index_end_t);
end
%===========================================================

%===========================================================
%Determine the satellite elevation angle.
%---------------------------------------------------------------------------------------------------------
la_dms =  [51 15 31.11582];
lo_dms = [-114 06 01.76988];
[XX_ref,Y_ref,Z_ref] = geodetic2ECEF(deg2rad(dms2degrees(la_dms)), ...
    deg2rad(dms2degrees(lo_dms)), ...
    1127.345, ...
    6378137, ...
    1/298.257223563); % WGS-84
x_ref = [XX_ref;Y_ref;Z_ref];
la_dms_true = [51 16 37.34162];
lo_dms_true = [-113 58 59.51154];
[X_base,Y_base,Z_base] = geodetic2ECEF(deg2rad(dms2degrees(la_dms_true)), ...
    deg2rad(dms2degrees(lo_dms_true)), ...
    1090.833, ...
    6378137, ...
    1/298.257223563); % WGS-84
x_base = [X_base,Y_base,Z_base];
ele = zeros(n_obs,1);
count = 1;
for i = 1:n_obs
    if(i <= count*12 && pr_s(i) ~= 0)
        ele(i) = elevation(x_ref',XS(i,:));
    elseif(i > count*12 && pr_s(i) ~= 0)
        count = count + 1;
        ele(i) = elevation(x_ref',XS(i,:));
    end
end
%===========================================================

% %===========================================================
% %Single Difference Iterative Least Square Formulation.
% %---------------------------------------------------------------------------------------------------------
% n_iteration = 0;
% XR0 = [0 0 0]; %initial value for reveiver position
% RES = cell(n_epoch,1);
% PRN = cell(n_epoch,1);
% n_sat_common = zeros(n_epoch,1);
% XR = zeros(n_epoch,3);
% dT = zeros(n_epoch,1);
% ELE = cell(n_epoch,1);
% EDOP = zeros(n_epoch,1);
% NDOP = zeros(n_epoch,1);
% VDOP = zeros(n_epoch,1);
% HDOP = zeros(n_epoch,1);
% PDOP = zeros(n_epoch,1);
% Cx_hat = zeros(n_epoch,3);
% 
% %Data Matching
% for i = 1:n_epoch
%     for k = 0:11
%         for j = 0:11
%             if (prn_s(12*i-11+k) == prn_rem(12*i-11+j))
%                 pr_s(12*i-11+k) = pr_rem(12*i-11+j);
%             end
%             if (prn_s(12*i-11+k) == prn_ref(12*i-11+j))
%                 pr_s_ref(12*i-11+k) = pr_ref(12*i-11+j);
%             end
%         end
%     end
% end
% 
% for i = 1:n_epoch
%     %find the common satellite PRN between remote and reference receiver
%     p1 = find(prn_rem((12*(i-1)+1):1:12*i));
%     for ii = 1:length(p1)
%         p1(ii) = p1(ii) + 12*(i-1);
%     end
%     p2 = find(prn_ref((12*(i-1)+1):1:12*i));
%     for iii = 1:length(p2)
%         p2(iii) = p2(iii) + 12*(i-1);
%     end
%     common_PRN = intersect(prn_rem(p1),...
%         prn_ref(p2));
%     PRN{i} = common_PRN;
%     n_sat_common(i) = length(PRN{i});
%     RES{i} = zeros(n_sat_common(i),1);
%     index = zeros(n_sat_common(i),1);
%     %index_ref = zeros(n_sat_common(i),1);
%     for j = 1:n_sat_common(i)
%         for k = 1:12
%             if prn_s(12*(i-1)+k) == PRN{i}(j)
%                 index(j) = 12*(i-1)+k;
%                 ELE{i}(j) = ele(12*(i-1)+k);
%             end
%         end
%     end
%     if(n_sat_common(i) >= min_sat)
%         while (n_iteration < 10)
%             [XR(i,:),dT(i),EDOP(i),NDOP(i),VDOP(i),HDOP(i),PDOP(i),Cx_hat(i,:),RES{i}] = ...
%                 SD_LS_adjustment(XR0,x_base,XS(index,:),pr_s(index),pr_s_ref(index));
% 
%             delta_x = XR0 - XR(i,:);
% 
%             %exit judgement
%             if (norm(delta_x(1:3)) < 1e-7)
%                 break;
%             end
% 
%             XR0 = XR(i,:);
%             n_iteration = n_iteration + 1;
%         end
%     end
%     n_iteration = 0;
% end
% %===========================================================

%===========================================================
%Double Difference Iterative Least Square Formulation.
%---------------------------------------------------------------------------------------------------------
n_iteration = 0;
XR0 = [0 0 0]; %initial value for reveiver position
RES = cell(n_epoch,1);
PRN = cell(n_epoch,1);
ELE = cell(n_epoch,1);
n_sat_common = zeros(n_epoch,1);
XR = zeros(n_epoch,3);
VDOP = zeros(n_epoch,1);
HDOP = zeros(n_epoch,1);
PDOP = zeros(n_epoch,1);
EDOP = zeros(n_epoch,1);
NDOP = zeros(n_epoch,1);
Cx_hat = zeros(n_epoch,3);

%Data Matching
for i = 1:n_epoch
    for k = 0:11
        for j = 0:11
            if (prn_s(12*i-11+k) == prn_rem(12*i-11+j))
                pr_s(12*i-11+k) = pr_rem(12*i-11+j);
            end
            if (prn_s(12*i-11+k) == prn_ref(12*i-11+j))
                pr_s_ref(12*i-11+k) = pr_ref(12*i-11+j);
            end
        end
    end
end

pivot = zeros(n_epoch,1);

for i = 1:n_epoch
    %find the common satellite PRN between remote and reference receiver
    p1 = find(prn_rem((12*(i-1)+1):1:12*i));
    for ii = 1:length(p1)
        p1(ii) = p1(ii) + 12*(i-1);
    end
    p2 = find(prn_ref((12*(i-1)+1):1:12*i));
    for iii = 1:length(p2)
        p2(iii) = p2(iii) + 12*(i-1);
    end
    common_PRN = intersect(prn_rem(p1),...
        prn_ref(p2));
    PRN{i} = common_PRN;
    n_sat_common(i) = length(PRN{i});
    RES{i} = zeros(n_sat_common(i),1);
    index = zeros(n_sat_common(i),1);
    
    for j = 1:n_sat_common(i)
        for k = 1:12
            if prn_s(12*(i-1)+k) == PRN{i}(j)
                index(j) = 12*(i-1)+k;
                ELE{i}(j) = ele(12*(i-1)+k);
            end
        end
    end
    
    max_ele = max(ele(index));
    pivot(i) = find(ele(index) == max_ele);
    
    if(n_sat_common(i) >= min_sat)
        while (n_iteration < 10)
            [XR(i,:),EDOP(i),NDOP(i),VDOP(i),HDOP(i),PDOP(i),Cx_hat(i,:),RES{i}] = ...
                DD_LS_adjustment_B(XR0,x_base,XS(index,:),pr_s(index),pr_s_ref(index),pivot(i));
            
            delta_x = XR0 - XR(i,:);
            
            %exit judgement
            if (norm(delta_x(1:3)) < 1e-7)
                break;
            end
            
            XR0 = XR(i,:);
            n_iteration = n_iteration + 1;
        end
    end
    n_iteration = 0;
end
%===========================================================

%===========================================================
%Compare with the Reference Solution.
%---------------------------------------------------------------------------------------------------------
diff_s = zeros(3,n_epoch);

for i = 1:n_epoch
    diff_s(:,i) = ECEF2ENU(x_ref,XR(i,:)');
end
%===========================================================

%===========================================================
%Accuracy Comparison Plot.
%---------------------------------------------------------------------------------------------------------
tt = 1:1:n_epoch;
ttt = 1:1:n_obs;
figure(1)
h1 = subplot(3,1,1);
plot(tt,diff_s(1,:),'-r',tt,ones(n_epoch,1)*mean(diff_s(1,:)),'-b',tt,Cx_hat(:,1),'-k',tt,-Cx_hat(:,1),'-k','LineWidth',1);
set(gca,'fontsize',14)
grid on
xlabel('Epoch [s]','FontSize', 14)
ylabel('East Error [m]','FontSize', 14)
title(h1,'East Error','FontSize', 14)
legend(h1,{'RMSE','Mean of RMSE','Estimated Accuracy'})
h2 = subplot(3,1,2);
plot(tt,diff_s(2,:),'-r',tt,ones(n_epoch,1)*mean(diff_s(2,:)),'-b',tt,Cx_hat(:,2),'-k',tt,-Cx_hat(:,2),'-k','LineWidth',1);
set(gca,'fontsize',14)
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
hold on
%===========================================================

%===========================================================
%Satellite DOP and total number plot.
%---------------------------------------------------------------------------------------------------------
figure(2)
plot(tt,EDOP,'b',tt,NDOP,'k',tt,VDOP,'r',tt,HDOP,'y',tt,PDOP,'g','LineWidth',2);
set(gca,'fontsize',14)
grid on
legend( 'EDOP','NDOP','VDOP','HDOP','PDOP');
xlabel('Epoch [s]','FontSize', 14)
title('DOP','FontSize', 14)

set(gcf,'units','points','position',[10,10,550,290])
hold on
%===========================================================

%===========================================================
%Residual Plot.
%---------------------------------------------------------------------------------------------------------
sat_type = intersect(prn_s,prn_rem);
sat_type = sat_type(sat_type ~= 0); %exclude 0
n_sat_type = length(sat_type);

figure(3)
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
p8 = p8(p8 ~= 0);t8 = t8(p8 ~= 0);e8 = e8(p8 ~= 0);
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
% plot(t7,p7,'d','Marker','.','Color',color_scheme(1,:),'markers',12);hold on;
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
set(gcf,'units','points','position',[10,10,660,320])
set(gca,'fontsize',14)
%===========================================================

%===========================================================
%Residual-Elevation Plot.
%---------------------------------------------------------------------------------------------------------
figure(4)
% plot(e7,p7,'d','Marker','.','Color',color_scheme(1,:),'markers',12);hold on;
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
ylim([-6 4])
title('Residuals of Observations - Satellite Elevation Angle','FontSize', 14)
set(gcf,'units','points','position',[10,10,660,320])
set(gca,'fontsize',14)
%===========================================================
