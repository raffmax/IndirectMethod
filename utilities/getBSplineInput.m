function [u,u_xi,u_T] = getBSplineInput(t,xi,T,nDeg,nSeg)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
nXi = numel(xi);

tau    = t/T;
nCtrl  = nSeg+nDeg; % number of control points
dKnot  = 1/nSeg;
% construct knot vector
knotVec = -nDeg*dKnot:dKnot:1+nDeg*dKnot;

u     = zeros(1,1);
u_tau = zeros(1,1);
u_xi = zeros(1,nXi);

xi_idx = 0; %index xi
for i = 0:nCtrl-1 
    % compute base spline
    B = baseSplineB(tau,knotVec,i,nDeg);
    % comput deriv. of base spline
    if nDeg>0
        BPrime = diffB(1,tau,knotVec,i,nDeg);
    end
    xi_idx  = xi_idx+1;
    xi_i = xi(xi_idx);
       
    u = u + xi_i*B;

    if nDeg>0
        u_tau = u_tau + xi_i*BPrime;
    end
    u_xi(1,xi_idx) = u_xi(1,xi_idx) + B;
end

% u_t  = u_tau/T;
u_T  = u_tau*(-t/(T^2));
end

function BPrime = diffB(order,t,knotVec,i,k)
    denom1 = knotVec(i+1+k)-knotVec(i+1);
    denom2 = knotVec(i+2+k)-knotVec(i+2);
    if denom1==0
        B_term1 = 0;
    else
        if order == 1
            B_term1 = k*baseSplineB(t,knotVec,i,k-1)/denom1;
        else
            B_term1 = k*diffB(order-1,t,knotVec,i,k-1)/denom1;
        end
    end
    if denom2==0
        B_term2 = 0;
    else
        if order == 1
            B_term2 = -k*baseSplineB(t,knotVec,i+1,k-1)/denom2;
        else
            B_term2 = -k*diffB(order-1,t,knotVec,i+1,k-1)/denom2;
        end
    end
   BPrime = B_term1+B_term2;
end

function B = baseSplineB(t,knotVec,i,k)
    % Cox–de Boor recursion 
    if k==0
        if t>=knotVec(i+1) && t<knotVec(i+2)
            B=1;
        else
            B=0;
        end
    else
        denom1 = knotVec(i+1+k)-knotVec(i+1);
        if denom1==0
            B_term1 = 0;
        else
            B_term1 = (t-knotVec(i+1))/denom1*baseSplineB(t,knotVec,i,k-1);
        end
        denom2 = knotVec(i+2+k)-knotVec(i+2);
        if denom2 == 0
            B_term2 = 0;
        else
            B_term2 = (knotVec(i+2+k)-t)/denom2*baseSplineB(t,knotVec,i+1,k-1);
        end
        B = B_term1+B_term2;
    end
end