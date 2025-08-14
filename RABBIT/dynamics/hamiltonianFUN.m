function H = hamiltonianFUN(x,p,q,u,gamma)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
% eval/call functions
f  = fFUN(x,u,gamma);
l  = lFUN(x,u); 

H = p'*f+q*l;
end