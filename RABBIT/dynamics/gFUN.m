function [g,g_x] = gFUN(x)

BTrafo = BTrafoAUTO(x);
W_sw   = W_swAUTO(x);
M      = MAUTO(x);

% Minimal coordinate discrete map
MCON = [M,-W_sw;-W_sw',zeros(2)];
res  = MCON\[M*BTrafo*[zeros(5),eye(5)]*x;0;0];
Delta__dq      = res(1:7);  % Discrete map for impact
Delta_flip__dq = blkdiag([0 1 0 0; 1 0 0 0; 0 0 0 1; 0 0 1 0],eye(3)) * Delta__dq;
Delta_min__dq  = [eye(5), zeros(5,2)] * Delta_flip__dq;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g = [[[0 1 0 0 0; 1 0 0 0 0; 0 0 0 1 0; 0 0 1 0 0; 0 0 0 0 1],zeros(5)]*x; Delta_min__dq];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dBTrafodx = dBTrafodxAUTO(x);
dW_swdx   = dW_swdxAUTO(x);
dMdx      = dMdxAUTO(x);

dresdx = zeros(9,10);
for i=1:10
    dMCONdxi = [dMdx(:,:,i),-dW_swdx(:,:,i);-dW_swdx(:,:,i)',zeros(2)];
    dxdxi    = zeros(10,1);
    dxdxi(i) = 1;
    dresdx(:,i) = -MCON\(dMCONdxi*res)...
                  +MCON\[dMdx(:,:,i)*BTrafo*[zeros(5),eye(5)]*x;0;0]...
                  +MCON\[M*dBTrafodx(:,:,i)*[zeros(5),eye(5)]*x;0;0]...
                  +MCON\[M*BTrafo*[zeros(5),eye(5)]*dxdxi;0;0];
end
dDelta__dqdx       = dresdx(1:7,:);
dxDelta_flip__dqdx = blkdiag([0 1 0 0; 1 0 0 0; 0 0 0 1; 0 0 1 0],eye(3)) * dDelta__dqdx;
dDelta_min__dqdx   = [eye(5), zeros(5,2)] * dxDelta_flip__dqdx;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g_x = [[[0 1 0 0 0; 1 0 0 0 0; 0 0 0 1 0; 0 0 1 0 0; 0 0 0 0 1],zeros(5)]; dDelta_min__dqdx];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end