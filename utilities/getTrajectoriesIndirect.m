function [t,x,u,cost,p] = getTrajectoriesIndirect(T,x0,p0,q,u0,gamma,vAVG,odeOpts)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

n = numel(x0);

[~,~,~,~,t,z] = flowIndirect(T,x0,p0,q,u0,odeOpts);
x = z(:,1:n);
yT = z(end,n+1)';
p = z(:,n+1+(1:n));
u = z(:,2*n+2:end);
cost = cFUN(T,x(end,:)',yT,gamma,vAVG);
end