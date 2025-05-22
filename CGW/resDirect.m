function [res,f] = resDirect(X,n,nXi,odeOpts)
% direct shooting implementation of r

% X  = [T,x0,xi,lambda,gamma,vAVG]
idx_T = 1;
idx_x0 = 1+(1:n);
idx_xi = 1+n+(1:nXi);
idx_lambda = 1+n+nXi+(1:n+2);
idx_gamma = 1+n+nXi+n+2+1;
idx_vAVG = 1+n+nXi+n+2+2;


T      = X(idx_T);
x0     = X(idx_x0);
xi     = X(idx_xi);
lambda = X(idx_lambda);
gamma  = X(idx_gamma); % slope
vAVG   = X(idx_vAVG); % avg speed

[zT,zT_T,zT_x0,zT_xi] = flowDirect(T,x0,xi,gamma,odeOpts);
xT = zT(1:n);
yT = zT(n+1);
xT_T  = zT_T(1:n,:);
yT_T  = zT_T(n+1,:);
xT_x0 = zT_x0(1:n,:);
yT_x0 = zT_x0(n+1,:);
xT_xi = zT_xi(1:n,:);
yT_xi = zT_xi(n+1,:);


[g,g_x] = gFUN(xT);
[hCon,hCon_T,hCon_x,~] = hFUN(T,xT,gamma,vAVG);
h = [g-x0;
     hCon];

Dh = [g_x*xT_T,g_x*xT_x0-eye(n),g_x*xT_xi;
      hCon_T+hCon_x*xT_T,hCon_x*xT_x0,hCon_x*xT_xi];

[c,c_T,c_x,c_y] = cFUN(T,xT,yT,gamma,vAVG);
f = c;
Df = [c_T+c_x*xT_T+c_y*yT_T,c_x*xT_x0+c_y*yT_x0,c_x*xT_xi+c_y*yT_xi];

res = [Df' + Dh'*lambda; h];

end

