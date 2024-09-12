function [t,x,u,cost] = getTrajectoriesDirect(T,x0,xi,gamma,vAVG,odeOpts)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

[zT,~,~,~,t,Z] = flowDirect(T,x0,xi,gamma,odeOpts);
yT = zT(end);

u = zeros(size(t))';
for i = 1:numel(t)
    if odeOpts.bezierVSbspline
        %[u,u_xi,u_T] = getPolyInput(t,xi);
        u(i) = getBSplineInput(t(i),xi,T,3,numel(xi)-3);
    else
        u(i) = getPolyInput(t(i),xi);
    end
end

x=Z(:,1:numel(x0));
cost = cFUN(T,x(end,:)',yT,gamma,vAVG);
end