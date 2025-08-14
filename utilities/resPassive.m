function [res,Dres] = resPassive(X,n,m,odeOpts)
% shooting implementation of r

% X  = [T,x0,xi,lambda,gamma,vAVG]
idx_T     = 1;
idx_x0    = 1+(1:n);
idx_gamma = 1+n+1;
idx_vAVG  = 1+n+2;

T     = X(idx_T);
x0    = X(idx_x0);
gamma = X(idx_gamma); % slope
vAVG  = X(idx_vAVG); % average speed

[xT,xT_T,xT_x0,xT_gamma] = flowPassive(T,x0,gamma,m,odeOpts);

[g,g_xT] = gFUN(xT);
[hCon,hCon_T,hCon_x,hCon_gamma] = hFUN(T,xT,gamma,vAVG);

res = [g-x0;
       hCon];

Dres = [g_xT*xT_T,g_xT*xT_x0-eye(n),g_xT*xT_gamma,zeros(n,1);...
        hCon_T+hCon_x*xT_T,hCon_x*xT_x0,hCon_gamma+hCon_x*xT_gamma,[0;-T]];


end

