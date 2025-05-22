function [res,f] = resIndirect(X,n,odeOpts)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

% X  = [T,x0,p0,lambda,gamma,vAVG]
idx_T = 1;
idx_x0 = 1+(1:n);
idx_p0 = 1+n+(1:n);
idx_q = 2+2*n;
idx_u0 = 2+2*n+1;
idx_lambda = 3+2*n+(1:2);
idx_gamma = 3+2*n+2+1;
idx_vAVG = 3+2*n+2+2;


T      = X(idx_T);
x0     = X(idx_x0);
p0     = X(idx_p0);
q     = X(idx_q);
u0     = X(idx_u0);
lambda = X(idx_lambda);
gamma  = X(idx_gamma); % slope
vAVG   = X(idx_vAVG); % avg speed

[xT,yT,pT,uT] = flowIndirect(T,x0,p0,q,u0,odeOpts);

[g,g_x] = gFUN(xT);
[h,h_T,h_x,~] = hFUN(T,xT,gamma,vAVG);


[c,c_T,c_x,c_y] = cFUN(T,xT,yT,gamma,vAVG);
f = c;

[~,~,f_u] = fFUN(x0,u0);
[~,~,l_u] = lFUN(x0,u0);
H0_u = f_u'*p0+l_u'*q;

HT = hamiltonianFUN(xT,pT,q,uT);

res = [g-x0;
       H0_u;
       q-c_y;
       h;
       HT+c_T+lambda'*h_T;
       g_x'*p0-pT+c_x'+h_x'*lambda];

end

