function [f,f_x,f_u] = fFUN(x,u,~)

f   = fAUTO(x,u);
f_x = f_xAUTO(x,u);
f_u = f_uAUTO(x,u);

end

