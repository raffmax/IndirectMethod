function [f,f_x,f_u,f_gamma] = fFUN(x,u,~)

f   = fAUTO(x,u);
f_x = f_xAUTO(x,u);
f_u = f_uAUTO(x,u);
f_gamma = zeros(4,1);

end

