function [zT,zT_T,zT_x0,zT_xi,t,Z] = flowDirect(T,x0,xi,gamma,odeOpts)
%UNTITLED Summary of this function goes here
%   z = [x;y]
%   Z = [z;z_T;z_x0;z_xi]

n   = numel(x0);
nXi = numel(xi);

odefun = @(t,Z) dZdtFUN(t,Z,n,xi,T,gamma,odeOpts.bezierVSbspline);

Z0 = [x0;0;...
      zeros(n,1);0;...
      reshape(eye(n),[n^2,1]);reshape(zeros(1,n),[n,1]);...
      reshape(zeros(n,nXi),[n*nXi,1]);reshape(zeros(1,nXi),[nXi,1])];
[t,Z] = ode15s(odefun,linspace(0,T,201),Z0,odeOpts);
zT = Z(end,1:(n+1))';
xT = zT(1:n);

if odeOpts.bezierVSbspline
    uT = getBSplineInput(T,xi,T,3,numel(xi)-3);
else
    uT = getPolyInput(T,xi);
end

zT_T  = [Z(end,n+1+(1:n))'+fFUN(xT,uT,gamma);Z(end,2*n+2)+lFUN(xT,uT)];
zT_x0 = [reshape(Z(end,2*n+2+(1:n^2))',[n,n]);...
         reshape(Z(end,2*n+2+n^2+(1:n)),[1,n])];
zT_xi = [reshape(Z(end,2*n+2+n^2+n+(1:n*nXi)),[n,nXi]);...
         reshape(Z(end,2*n+2+n^2+n+n*nXi+(1:nXi)),[1,nXi])];
end

function dZdt = dZdtFUN(t,Z,n,xi,T,gamma,bezierVSbspline)
nXi = numel(xi);

x = Z(1:n); % dynamic states
%y = Z(n+1); % integrand
x_T  = Z(n+1+(1:n));
% y_T  = Z(2*n+2);
x_x0 = reshape(Z(2*n+2+(1:n^2)),[n,n]);
%y_x0 = reshape(Z(2*n+2+n^2+(1:n)),[1,n]);
x_xi = reshape(Z(2*n+2+n^2+n+(1:n*nXi)),[n,nXi]);
%y_xi = reshape(Z(2*n+2+n^2+n+n*nXi+(1:nXi)),[1,nXi]);

if bezierVSbspline
    [u,u_xi,u_T] = getBSplineInput(t,xi,T,3,numel(xi)-3);
else
    [u,u_xi,u_T] = getPolyInput(t,xi);
end

[f,f_x,f_u] = fFUN(x,u,gamma);
[l,l_x,l_u] = lFUN(x,u);

dx_x0dt = f_x*x_x0;
dy_x0dt = l_x*x_x0;

dx_Tdt = f_x*x_T+f_u*u_T;
dy_Tdt = l_x*x_T+l_u*u_T;

f_xi = f_u*u_xi;
l_xi = l_u*u_xi;
dx_xidt = f_x*x_xi+f_xi;
dy_xidt = l_x*x_xi+l_xi;

dZdt = [f;l;...
      dx_Tdt;dy_Tdt;...
      reshape(dx_x0dt,[n^2,1]);reshape(dy_x0dt,[n,1]);...
      reshape(dx_xidt,[n*nXi,1]);reshape(dy_xidt,[nXi,1])];
end