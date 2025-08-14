function [bifInd,Dres_min] = bifurcationIndicator(T,rFun,x_equi)
    [~,Dres] = rFun([T;x_equi;0;0]);
    % Laypunov-Schmidt reduction: Dres_min
    % remove event e=0 and period T
    Dres_min = Dres([1:end-2,end],2:end);
    % further isolating/augmenting the Jacobian by tangent vector,
    % for which: q5 = -gamma
    Dres_min = [Dres_min; zeros(1,4) 1 zeros(1,5) 1 0];
    bifInd = det(Dres_min);
end