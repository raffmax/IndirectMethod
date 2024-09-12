function [l,l_x,l_u,l_ux] = lFUN(x,u)
% stage cost
l = 0.5*u^2;
l_x = [0,0,0,0];
l_u = u;

l_ux = [0,0,0,0];

end
