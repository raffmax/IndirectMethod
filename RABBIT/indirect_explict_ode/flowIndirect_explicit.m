function [xT,yT,pT,uT,t,Z] = flowIndirect_explicit(T,x0,p0,q,gamma,odeOpts)
%UNTITLED Summary of this function goes here
%   z = [x;y;p]

n   = numel(x0);
N   = 201;

odefun = @(t,Z) dZdtFUN(t,Z,n,q,gamma);

Z0 = [x0;0;p0];
[t,Z] = odeOpts.SolverType(odefun,linspace(0,T,N),Z0,odeOpts);
xT = Z(end,1:n)';
yT = Z(end,n+1)';
pT = Z(end,n+1+(1:n))';

B_min = [eye(4);0 0 0 0];
M_min = M_minAUTO(xT);
f_u   = [zeros(5,4);M_min\B_min];
uT     = -f_u'*pT/q;

m  = numel(uT);
x_ = Z(:,1:n);
y_ = Z(:,n+1);
p_ = Z(:,n+1+(1:n));
Z  = zeros(N,2*n+1+m);
for i=1:N
    x = x_(i,:)';
    y = y_(i,:)';
    p = p_(i,:)';
    B_min = [eye(4);0 0 0 0];
    M_min = M_minAUTO(x);
    f_u   = [zeros(5,4);M_min\B_min];
    u     = -f_u'*p/q;
    Z(i,:)= [x',y',p',u'];
end
end

function dZdt = dZdtFUN(~,Z,n,q,gamma)
x = Z(1:n); % dynamic states
%y = Z(n+1); % integrand
p = Z(n+1+(1:n)); % co-states

% B_min = [eye(4);0 0 0 0];
% M_min = M_minAUTO(x);
% f_u   = [zeros(5,4);M_min\B_min];
% u     = -f_u'*p/q;
% 0 = f_u'*p+l_u'*q;
B_min = [eye(4);0 0 0 0];
M_min = M_minAUTO(x);
f_u   = [zeros(5,4);M_min\B_min];
u     = -f_u'*p/q;

[f,f_x,~] = fFUN(x,u,gamma);
[l,l_x,~] = lFUN(x,u);

dxdt = f;
dydt = l;
dpdt = -f_x'*p-l_x'*q;

dZdt = [dxdt;dydt;dpdt];
end