function [l,l_x,l_u,l_ux] = lFUN(x,u)
% stage cost
l = 0.5*u'*u;
l_x = zeros(1,10);
l_u = u';

l_ux = zeros(4,10);


end
