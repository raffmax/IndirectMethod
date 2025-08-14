function [zT,zT_T,zT_x0,zT_xi,t,Z] = flowDirect(T,x0,xi,gamma,m,odeOpts)
%UNTITLED Summary of this function goes here
%   z = [x;y]
%   Z = [z;z_T;z_x0;z_xi]
%   m = dim(u)

n   = numel(x0);
nXi = numel(xi)/m;

odefun = @(t,Z) dZdtFUN(t,Z,n,m,xi,T,gamma,odeOpts.bezierVSbspline);

Z0 = [x0;0;...
      zeros(n,1);0;...
      reshape(eye(n),[n^2,1]);reshape(zeros(1,n),[n,1]);...
      reshape(zeros(n,nXi*m),[n*nXi*m,1]);reshape(zeros(1,nXi*m),[nXi*m,1])];
[t,Z] = odeOpts.SolverType(odefun,linspace(0,T,201),Z0,odeOpts);
zT = Z(end,1:(n+1))';
xT = zT(1:n);

uT    = zeros(m,1);
for i = 1:m
    if odeOpts.bezierVSbspline
        ui = getBSplineInput(T,xi(nXi*(i-1)+(1:nXi)),T,3,nXi-3);
    else
        ui = getPolyInput(T,xi(nXi*(i-1)+(1:nXi)));
    end
    uT(i) = ui;
end


zT_T  = [Z(end,n+1+(1:n))'+fFUN(xT,uT,gamma);Z(end,2*n+2)+lFUN(xT,uT)];
zT_x0 = [reshape(Z(end,2*n+2+(1:n^2))',[n,n]);...
         reshape(Z(end,2*n+2+n^2+(1:n)),[1,n])];
zT_xi = [reshape(Z(end,2*n+2+n^2+n+(1:n*nXi*m)),[n,nXi*m]);...
         reshape(Z(end,2*n+2+n^2+n+n*nXi*m+(1:nXi*m)),[1,nXi*m])];
end

function dZdt = dZdtFUN(t,Z,n,m,xi,T,gamma,bezierVSbspline)
nXi = numel(xi)/m;

x = Z(1:n); % dynamic states
%y = Z(n+1); % integrand
x_T  = Z(n+1+(1:n));
% y_T  = Z(2*n+2);
x_x0 = reshape(Z(2*n+2+(1:n^2)),[n,n]);
%y_x0 = reshape(Z(2*n+2+n^2+(1:n)),[1,n]);
x_xi = reshape(Z(2*n+2+n^2+n+(1:n*nXi*m)),[n,nXi*m]);
%y_xi = reshape(Z(2*n+2+n^2+n+n*nXi+(1:nXi)),[1,nXi]);

u    = zeros(m,1);
u_xi = zeros(m,nXi*m);
u_T  = zeros(m,1);
for i = 1:m
    if bezierVSbspline
        [ui,ui_xi,ui_T] = getBSplineInput(t,xi(nXi*(i-1)+(1:nXi)),T,3,nXi-3);
    else
        [ui,ui_xi,ui_T] = getPolyInput(t,xi(nXi*(i-1)+(1:nXi)));
    end
    u(i)      = ui;
    u_xi(i,:) = [zeros(1,nXi*(i-1)),ui_xi,zeros(1,nXi*(m-i))];
    u_T(i)    = ui_T;
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
      reshape(dx_xidt,[n*nXi*m,1]);reshape(dy_xidt,[nXi*m,1])];
end