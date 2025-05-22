function [c,c_T,c_x,c_y] = cFUN(T,xT,yT,gamma,vAVG)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

c = yT/(T*vAVG);
c_T = -yT/(T^2*vAVG);
c_x = [0,0,0,0];
c_y = 1/(T*vAVG);
end

