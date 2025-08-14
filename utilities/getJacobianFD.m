function [f,jacobian] = getJacobianFD(fcn,x,stepSize)
%GETJACOBIANFD Finite Differences on given function
%   using central finite differences

% get f(x)
f  = fcn(x);
% get dimensions
nF = length(f);
nX = length(x);

jacobian = zeros(nF,nX);
h        = stepSize;
for iCol = 1:nX
       iX       = zeros(nX,1);
       iX(iCol) = 1;
       jacobian(:,iCol) = (fcn(x+h*iX)-fcn(x-h*iX))/(2*h);
end

end

