%Read observation and satellite files
[prn_s,t_s,x,y,z,xv,yv,zv] = readsat('Satellites.sat');
[n_obs,prn_o,t_o,pr,cp_l1,doppler_l1,cp_l2] = readobs('RemoteL1L2.obs');
XS = [x y z]; %Satellites Corrdinate (ECEF)

n_epoch = n_obs / 12; %total number of epoches
n_sat = zeros(n_epoch,1); %the number of satelltes on each epoch

PRN = [7 8 9 11 15 17 18 19 22 24 26 27 28];
AZ = zeros(n_epoch,13);
EL = zeros(n_epoch,13);

la_dms =  [51 15 31.11582];
lo_dms = [-114 06 01.76988];
[XX_ref,Y_ref,Z_ref] = geodetic2ECEF(deg2rad(dms2degrees(la_dms)), ...
    deg2rad(dms2degrees(lo_dms)), ...
    1127.345, ...
    6378137, ...
    1/298.257223563); % WGS-84
x_ref = [XX_ref;Y_ref;Z_ref];

count = 1;
ele = zeros(n_obs,1);
az = zeros(n_obs,1);
for i = 1:n_obs
    if(i <= count*12)
        [az(i),ele(i)] = elee(x_ref',XS(i,:));
    elseif(i > count*12)
        count = count + 1;
        [az(i),ele(i)] = elee(x_ref',XS(i,:));
    end
end

count_7 = 1;count_8 = 1;count_9 = 1;count_11 = 1;count_15 = 1;count_17 = 1;
count_18 = 1;count_19 = 1;count_22 = 1;count_24 = 1;count_26 = 1;
count_27 = 1;count_28 = 1;


for i = 1:n_obs
    if prn_s(i)  == 7
        AZ(count_7,1) = az(i);
        EL(count_7,1) = ele(i);
        count_7 = count_7 + 1;
    elseif prn_s(i) == 8
        AZ(count_8,2) = az(i);
        EL(count_8,2) = ele(i);
        count_8 = count_8 + 1;
    elseif prn_s(i) == 9
        AZ(count_9,3) = az(i);
        EL(count_9,3) = ele(i);
        count_9 = count_9 + 1;
    elseif prn_s(i) == 11
        AZ(count_11,4) = az(i);
        EL(count_11,4) = ele(i);
        count_11 = count_11 + 1;
    elseif prn_s(i) == 15
        AZ(count_15,5) = az(i);
        EL(count_15,5) = ele(i);
        count_15 = count_15 + 1;
    elseif prn_s(i) == 17
        AZ(count_17,6) = az(i);
        EL(count_17,6) = ele(i);
        count_17 = count_17 + 1;
    elseif prn_s(i) == 18
        AZ(count_18,7) = az(i);
        EL(count_18,7) = ele(i);
        count_18 = count_18 + 1;
    elseif prn_s(i) == 19
        AZ(count_19,8) = az(i);
        EL(count_19,8) = ele(i);
        count_19 = count_19 + 1;
    elseif prn_s(i) == 22
        AZ(count_22,9) = az(i);
        EL(count_22,9) = ele(i);
        count_22 = count_22 + 1;
    elseif prn_s(i) == 24
        AZ(count_24,10) = az(i);
        EL(count_24,10) = ele(i);
        count_24 = count_24 + 1;
    elseif prn_s(i) == 26
        AZ(count_26,11) = az(i);
        EL(count_26,11) = ele(i);
        count_26 = count_26 + 1;
    elseif prn_s(i) == 27
        AZ(count_27,12) = az(i);
        EL(count_27,12) = ele(i);
        count_27 = count_27 + 1;
    elseif prn_s(i) == 28
        AZ(count_28,13) = az(i);
        EL(count_28,13) = ele(i);
        count_28 = count_28 + 1;
    end
end

for i = 1:n_epoch
    if EL(i) == 0
        EL(i) = 90;
    end
end

h = skyPlot(AZ',EL',PRN);


function hpol = skyPlot(varargin)
%Function plots "sky view" from the receiver perspective.
%
%h = skyPlot(AZ, EL, PRN, line_style)
%
%   Inputs:
%       AZ              - contains satellite azimuth angles. It is a 2D
%                       matrix. One line contains data of one satellite.
%                       The columns are the calculated azimuth values.
%       EL              - contains satellite elevation angles. It is a 2D
%                       matrix. One line contains data of one satellite.
%                       The columns are the calculated elevations.
%       PRN             - a row vector containing PRN numbers of the
%                       satellites.
%       line_style      - line style of the plot. The same style will be
%                       used to plot all satellite positions (including
%                       color).
%   Outputs:
%       h               - handle to the plot

%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
%
% Copyright (C) Darius Plausinaitis and Kristin Larson
% Written by Darius Plausinaitis and Kristin Larson
%--------------------------------------------------------------------------
%This program is free software; you can redistribute it and/or
%modify it under the terms of the GNU General Public License
%as published by the Free Software Foundation; either version 2
%of the License, or (at your option) any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License
%along with this program; if not, write to the Free Software
%Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
%USA.
%--------------------------------------------------------------------------

%CVS record:
%$Id: skyPlot.m,v 1.1.2.5 2006/08/18 11:41:57 dpl Exp $

%% Check arguments and sort them ==========================================
[hAxis, args, nargs] = axescheck(varargin{:});

if nargs < 3 || nargs > 4
    error('Requires 3 or 4 data arguments.')
elseif nargs == 3
    [az, el, prn]   = deal(args{1:3});
    line_style      = 'auto';
else
    [az, el, prn, line_style] = deal(args{1:4});
end

if ischar(az) || ischar(el) || ischar(prn)
    error('AZ and EL must be numeric.');
end

if ~isequal(size(az), size(el))
    error('AZ and EL must be same size.');
end

%% Prepare axis ===========================================================
hAxis = newplot(hAxis);

%--- Get x-axis text color so grid is in same color -----------------------
tc = get(hAxis, 'xcolor');

hold(hAxis, 'on');

%--- Plot white background ------------------------------------------------
rectangle('position', [-90, -90, 180, 180], ...
    'Curvature', [1 1], ...
    'facecolor', 'white', ...
    'edgecolor', tc);

%% Plot spokes ============================================================

%--- Find spoke angles ----------------------------------------------------
% Only 6 lines are needed to divide circle into 12 parts
th = (1:6) * 2*pi / 12;

%--- Convert spoke end point coordinate to Cartesian system ---------------
cst = cos(th); snt = sin(th);
cs = [cst; -cst];
sn = [snt; -snt];

%--- Plot the spoke lines -------------------------------------------------
line(90*sn, 90*cs, 'linestyle', ':', 'color', tc, 'linewidth', 0.5, ...
    'handlevisibility', 'off');

%% Annotate spokes in degrees =============================================
rt = 1.1 * 90;

for i = 1:max(size(th))
    
    %--- Write text in the first half of the plot -------------------------
    text(rt*snt(i), rt*cst(i), int2str(i*30), ...
        'horizontalalignment', 'center', 'handlevisibility', 'off');
    
    if i == max(size(th))
        loc = int2str(0);
    else
        loc = int2str(180 + i*30);
    end
    
    %--- Write text in the opposite half of the plot ----------------------
    text(-rt*snt(i), -rt*cst(i), loc, ...
        'handlevisibility', 'off', 'horizontalalignment', 'center');
end

%% Plot elevation grid ====================================================

%--- Define a "unit" radius circle ----------------------------------------
th = 0 : pi/50 : 2*pi;
xunit = cos(th);
yunit = sin(th);

%--- Plot elevation grid lines and tick text ------------------------------
for elevation = 0 : 15 : 90
    elevationSpherical = 90*cos((pi/180) * elevation);
    
    line(yunit * elevationSpherical, xunit * elevationSpherical, ...
        'lineStyle', ':', 'color', tc, 'linewidth', 0.5, ...
        'handlevisibility', 'off');
    
    text(0, elevationSpherical, num2str(elevation), ...
        'BackgroundColor', 'white', 'horizontalalignment','center', ...
        'handlevisibility', 'off');
end

%--- Set view to 2-D ------------------------------------------------------
view(0, 90);

%--- Set axis limits ------------------------------------------------------
%save some space for the title
axis([-95 95 -90 101]);

%% Transform elevation angle to a distance to the center of the plot ------
elSpherical = 90*cos(el * pi/180);

%--- Transform data to Cartesian coordinates ------------------------------
yy = elSpherical .* cos(az * pi/180);
xx = elSpherical .* sin(az * pi/180);

%% Plot data on top of the grid ===========================================

if strcmp(line_style, 'auto')
    %--- Plot with "default" line style -----------------------------------
    hpol = plot(hAxis, xx', yy', '.-');
else
    %--- Plot with user specified line style ------------------------------
    % The same line style and color will be used for all satellites
    hpol = plot(hAxis, xx', yy', line_style);
end

%--- Mark the last position of the satellite ------------------------------
plot(hAxis, xx(:,1)', yy(:,1)', 'o', 'MarkerSize', 10);

%--- Place satellite PRN numbers at the latest position -------------------
for i = 1:length(prn)
    if(prn(i) ~= 0)
        % The empthy space is used to place the text a side of the last
        % point. This solution results in constant offset even if a zoom
        % is used.
        text(xx(i, 1), yy(i, 1), ['  ', int2str(prn(i))], 'color', 'b');
    end
end

%--- Make sure both axis have the same data aspect ratio ------------------
axis(hAxis, 'equal');

%--- Switch off the standard Cartesian axis -------------------------------
axis(hAxis, 'off');
end

function [Az,El] = elee(XR, XS)
[phi, lam] = ECEF2geodetic(XR(1), XR(2), XR(3));

%rotation matrix from global to local reference system
R = [-sin(lam) cos(lam) 0;
    -sin(phi)*cos(lam) -sin(phi)*sin(lam) cos(phi);
    +cos(phi)*cos(lam) +cos(phi)*sin(lam) sin(phi)];

local = R*(XS-XR)';
E = local(1);
N = local(2);
U = local(3);
dis = sqrt(E^2+N^2);

if dis < 1.e-20
    Az = 0;
    El = 90;
else
    Az = rad2deg(atan2(E,N));
    El = rad2deg(atan2(U,dis));
end

if Az < 0
    Az = Az + 360;
end
end






















