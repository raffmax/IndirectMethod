function [xT,yT,pT,uT,t,Z] = flowIndirect(T,x0,p0,q,u0,odeOpts)
%UNTITLED Summary of this function goes here
%   z = [x;y;p]

n   = numel(x0);

% add mass matrix to solve DAE of index 1
odeOpts.Mass = blkdiag(eye(2*n+1),0);

odefun = @(t,Z) dZdtFUN(t,Z,n,q);

Z0 = [x0;0;p0;u0];
[t,Z] = ode15s(odefun,linspace(0,T,201),Z0,odeOpts);
xT = Z(end,1:n)';
yT = Z(end,n+1)';
pT = Z(end,n+1+(1:n))';
uT = Z(end,2*n+2:end)';
end

function dZdt = dZdtFUN(~,Z,n,q)
x = Z(1:n); % dynamic states
%y = Z(n+1); % integrand
p = Z(n+1+(1:n)); % co-states

u = Z(2*n+2:end);

[f,f_x,f_u] = fFUN(x,u);
[l,l_x,l_u] = lFUN(x,u);

dxdt = f;
dydt = l;
dpdt = -f_x'*p-l_x'*q;
dudt = f_u'*p+l_u'*q;

dZdt = [dxdt;dydt;dpdt;dudt];
end