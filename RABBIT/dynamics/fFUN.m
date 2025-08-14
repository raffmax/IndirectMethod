function [f,f_x,f_u,f_gamma] = fFUN(x,u,gamma)

B_min = [eye(4);0 0 0 0];
M_min = M_minAUTO(x);
c_min = c_minAUTO(x,gamma);
G_min = G_minAUTO(x,gamma);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = [x(6:10); M_min \ (B_min * u - G_min - c_min)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dM_mindx = dM_mindxAUTO(x);
dc_mindx = dc_mindxAUTO(x,gamma);
dG_mindx = dG_mindxAUTO(x,gamma);
dc_mindgamma = dc_mindgammaAUTO(x,gamma);
dG_mindgamma = dG_mindgammaAUTO(x,gamma);

f_dq = f(6:10);
dM_mindx__f_dq = zeros(5,10);
for i=1:10
    dM_mindx__f_dq(:,i) = dM_mindx(:,:,i)*f_dq;
end
df_qdx  = [zeros(5),eye(5)];
df_dqdx = M_min\(-dM_mindx__f_dq - dG_mindx - dc_mindx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_x     = [df_qdx;df_dqdx];
f_u     = [zeros(5,4);M_min\B_min];
f_gamma = [zeros(5,1);M_min\(-dG_mindgamma-dc_mindgamma)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

