%Read Satellite Binary File, with each element is 64-bit double
%Each second of data contains 12 channels. Each channel can hold coordinates for one satellite.
%Channels without satellites will contain zeroes. The format for one channel is as follows: 
%                                                       <Satellite PRN>
%                                           <GPS time of observation (seconds)> 
%                                               <X Satellite Coordinate (m)> 
%                                               <Y Satellite Coordinate (m)> 
%                                               <Z Satellite Coordinate (m)> 
%                                               <X Satellite Velocity (m/s)> 
%                                               <Y Satellite Velocity (m/s)>
%                                               <Z Satellite Velocity (m/s)> 
%Script Written by Cheng.

function [prn,t_0,x,y,z,xv,yv,zv] = readsat(file)
fileID = fopen(file,'r');
onebyte = fread(fileID,100000000,'double');

%disp(onebyte);

n = numel(onebyte); %number of elements in the file, can be divided by 8
m = n/8;

prn = zeros(m,1,'double');
t_0 = zeros(m,1,'double');
x = zeros(m,1,'double');
y = zeros(m,1,'double');
z = zeros(m,1,'double');
xv = zeros(m,1,'double');
yv = zeros(m,1,'double');
zv = zeros(m,1,'double');

count = 1;

XX = [0 ;0 ;0];

for i = 1:m
    if count < n
        prn(i) = onebyte(count);    count = count +1;
        t_0(i) = onebyte(count);    count = count +1;
        x(i) = onebyte(count);     count = count +1;
        y(i) = onebyte(count);      count = count +1;
        z(i) = onebyte(count);     count = count +1;
        xv(i) = onebyte(count);      count = count +1;
        yv(i) = onebyte(count);    count = count +1;
        zv(i) = onebyte(count);    count = count +1;
    end
end

end