function [u,u_xi,u_T] = getPolyInput(t,xi)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
nXi = numel(xi);
nT  = numel(t);

basis = ones(nT,nXi);

% % monom basis
% nM = nXi-1;
% for iM=0:nM
%     basis(:,iM+1) = t(:).^iM;
% end

% bernsein basis
nB = nXi-1;
for iB=0:nB
    basis(:,iB+1) = nchoosek(nB,iB).*(t(:).^iB).*((1-t(:)).^(nB-iB));
end

% % hermite basis
% nH = nXi-1;
% if nH>0
%     basis(:,2) = 2*t(:);
% end
% if nH>1
%     for iH=2:nH
%         basis(:,iH+1) = 2*t(:).*basis(:,iH)-2*(iH-1)*basis(:,iH-1);
%     end
% end

u = basis*xi;
u_xi = basis;
u_T = 0;
end