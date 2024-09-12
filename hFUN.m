function [h,h_T,h_x,h_gamma] = hFUN(T,xT,gamma,vAVG)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
h = [xT(1)+xT(2)+2*gamma;2*sin(xT(1)+gamma)-vAVG*T];
h_T = [0;-vAVG];
h_x = [1 1 0 0;2*cos(xT(1)+gamma) 0 0 0];
h_gamma = [2;2*cos(xT(1)+gamma)];
end