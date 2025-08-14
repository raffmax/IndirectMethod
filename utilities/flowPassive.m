function [xT,xT_T,xT_x0,xT_gamma,t,X] = flowPassive(T,x0,gamma,m,odeOpts)
%UNTITLED Summary of this function goes here
%   X = [x;x_x0;x_gamma]
%   m = dim(u)

n   = numel(x0);

odefun = @(t,X) dXdtFUN(t,X,n,m,gamma);

X0 = [x0;...
      reshape(eye(n),[n^2,1]);...
      zeros(n,1)];
[t,X] = odeOpts.SolverType(odefun,linspace(0,T,201),X0,odeOpts);
xT = X(end,1:n)';

xT_T     = fFUN(xT,zeros(m,1),gamma);
xT_x0    = reshape(X(end,n+(1:n^2))',[n,n]);
xT_gamma = X(end,n+n^2+(1:n))';
end

function dXdt = dXdtFUN(~,X,n,m,gamma)

x       = X(1:n); % dynamic states
x_x0    = reshape(X(n+(1:n^2)),[n,n]);
x_gamma = X(n+n^2+(1:n));

[f,f_x,~,f_gamma] = fFUN(x,zeros(m,1),gamma);

dx_x0dt = f_x*x_x0;

dx_gammadt = f_x*x_gamma+f_gamma;

dXdt = [f;...
        reshape(dx_x0dt,[n^2,1]);...
        dx_gammadt];
end