function contData = indirectContinuation(z_init,rFun,sigma_end,n,fDstepSize,rootFunctionTolerance,h)
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here

opts = [];
opts.Grad1 = false;
opts.aimOnTarget = true;
opts.idxConPar = numel(z_init);
opts.FiniteDifferenceStepSize = fDstepSize;
opts.MaxIterations = 20;
opts.FunctionTolerance = rootFunctionTolerance;
opts.StepTolerance = 1e-11;
[z,fval,~,output,jac] = NewtonsMethod(rFun,z_init,opts);

output.t  = getTangent(jac);
simJacobian = jac(:,1:opts.idxConPar-1);
augJacobian = [jac;output.t'];
contData.simJacobian = zeros(size(simJacobian,1),size(simJacobian,2),1e4);
contData.augJacobian = zeros(size(augJacobian,1),size(augJacobian,2),1e4);
contData.simJacDet = det(simJacobian);
contData.augJacDet = det(augJacobian);
contData.simJacobian(:,:,1) = simJacobian;
contData.augJacobian(:,:,1) = augJacobian;
contData.T = z(1);
contData.x0 = z(1+(1:n));
contData.p0 = z(1+n+(1:n));
contData.q = z(2*n+2);
contData.u0 = z(2*n+2+1);
contData.lambda = z(3+2*n+(1:2));
contData.sigma = z(end);
[~,cost] = rFun(z);
contData.cost = cost;


opts.aimOnTarget = false; 
endLoop = false;
kLoop = 1;
direction = sign(sigma_end-z_init(opts.idxConPar));
d = direction*sign(output.t(opts.idxConPar)); %direction
while ~endLoop
    kLoop = kLoop+1;
    % predictor step
    z_new = z+h*d*output.t; 
    % corrector step
    [z_new,fval,exitflag,output,jac] = NewtonsMethod(rFun,z_new,opts);
    if exitflag==-1
        break
    end
    if direction*z_new(end)>direction*sigma_end
        opts.aimOnTarget = true;
        endLoop = true;
        z_new = z_new-(z_new(end)-sigma_end)*output.t/output.t(end);
        [z,fval,~,output,jac] = NewtonsMethod(rFun,z_new,opts);
    else
        z = z_new;
    end
    simJacobian = jac(:,1:opts.idxConPar-1);
    augJacobian = [jac;output.t'];
    contData.simJacDet = [contData.simJacDet,det(simJacobian)];
    contData.augJacDet = [contData.augJacDet,det(augJacobian)];
    contData.simJacobian(:,:,kLoop) = simJacobian;
    contData.augJacobian(:,:,kLoop) = augJacobian;
    contData.T = [contData.T,z(1)];
    contData.x0 = [contData.x0,z(1+(1:n))];
    contData.p0 = [contData.p0,z(1+n+(1:n))];
    contData.q = [contData.q,z(2+2*n)];
    contData.u0 = [contData.u0,z(2*n+2+1)];
    contData.lambda = [contData.lambda,z(3+2*n+(1:2))];
    contData.sigma = [contData.sigma,z(end)];
    [~,cost] = rFun(z);
    contData.cost = [contData.cost,cost];
end
contData.simJacobian = contData.simJacobian(:,:,1:kLoop);
contData.augJacobian = contData.augJacobian(:,:,1:kLoop);
end