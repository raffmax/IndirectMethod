function [t,x] = getTrajectoriesPassive(T,x0,gamma,m,odeOpts)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

n = numel(x0);

[~,~,~,~,t,z] = flowPassive(T,x0,gamma,m,odeOpts);
x = z(:,1:n);
end