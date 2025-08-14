function [h,h_T,h_x,h_gamma] = hFUN(T,xT,~,vAVG)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
pos_swingfoot   = posSwingFootAUTO(xT);
pos_swingfoot_x = dposSwingFootdxAUTO(xT);
h = [pos_swingfoot(2);pos_swingfoot(1)-vAVG*T];
h_T = [0;-vAVG];
h_x = [pos_swingfoot_x(2,:);pos_swingfoot_x(1,:)];
h_gamma = [0;0];
end