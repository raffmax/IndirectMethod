function [t,x,u,cost] = getTrajectoriesDirect(T,x0,xi,gamma,m,vAVG,odeOpts)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

[zT,~,~,~,t,Z] = flowDirect(T,x0,xi,gamma,m,odeOpts);
yT = zT(end);

Xi = reshape(xi,m,numel(xi)/m);
u  = zeros(m,numel(t));
for j = 1:m
for i = 1:numel(t)
    if odeOpts.bezierVSbspline
        %[u,u_xi,u_T] = getPolyInput(t,xi);
        u(j,i) = getBSplineInput(t(i),Xi(:,j),T,3,numel(Xi(:,j))-3);
    else
        u(j,i) = getPolyInput(t(i),Xi(:,j));
    end
end
end

x=Z(:,1:numel(x0));
cost = cFUN(T,x(end,:)',yT,gamma,vAVG);
end