%Read Observation Binary File, with each element is 64-bit double
%Each second of data contains 12 channels. Each channel can hold observations for one satellite.
%Channels without observations will contain zeroes. The format for one channel is as follows:
%                                              <Satellite PRN>
%                                    <GPS time of observation (seconds)>
%                                           <C/A code pseudorange (m)>
%                                 <L1 carrier phase measurement (L1 cycles)>
%                                         <Doppler of L1 carrier (Hertz)>
%                                 <L2 carrier phase measurement (L2 cycles)>
%Script Written by Cheng Huang.

function [m,prn,t_0,prange,cp_l1,doppler_l1,cp_l2] = readobs(file)
fileID = fopen(file,'r');
onebyte = fread(fileID,10000000,'double');

%disp(onebyte);

n = numel(onebyte); %number of elements in the file, can be divided by 6
m = n/6;
prn = zeros(m,1,'double');
t_0 = zeros(m,1,'double');
prange = zeros(m,1,'double');
cp_l1 = zeros(m,1,'double');
doppler_l1 = zeros(m,1,'double');
cp_l2 = zeros(m,1,'double');

count = 1;

for i = 1:m
    if count < n
        prn(i) = onebyte(count);    count = count +1;
        t_0(i) = onebyte(count);    count = count +1;
        prange(i) = onebyte(count);     count = count +1;
        cp_l1(i) = onebyte(count);      count = count +1;
        doppler_l1(i) = onebyte(count);     count = count +1;
        cp_l2(i) = onebyte(count);      count = count +1;
    end
end

end


