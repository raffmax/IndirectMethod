%% set up
addpath('utilities');
addpath('dynamics');
addpath('animations');

n   = 4; % state dimension
nXi = 4; % number of control parameters

% passive solution at vAVG=0.1
load passiveGait_01.mat
vAVG = passiveGait_01.vAVG;

gamma_end = 0;

fDstepSize = 1e-9; % finite difference step size
rootFunctionTolerance = 1e-8;
odeOpts = odeset('RelTol',1e-9,'AbsTol',1e-10); % ode solver
odeOpts.bezierVSbspline = 1; % Bezier=0, B-Spline=1
h = 1e-2; % fixed step-size in continuation

% X = [T,x0,xi,lambda,gamma,vAVG]
% X = [z;vAVG]
X0 = [passiveGait_01.T;
      passiveGait_01.x0;
      zeros(nXi,1);
      zeros(n+2,1);
      passiveGait_01.gamma;
      vAVG];

rFun = @(z) resDirect([z;vAVG],n,nXi,odeOpts); % fix vAVG
% z = [T,x0,xi,lambda,gamma]
z_init = X0(1:end-1);

contData = directContinuation(z_init,rFun,gamma_end,n,nXi,...
                              fDstepSize,rootFunctionTolerance,h);

contData.gamma = contData.sigma;
contData.vAVG  = vAVG*ones(size(contData.T));
contData = rmfield(contData, 'sigma');


%% plots
idx = [1,find(diff(contData.strictMin<0)),numel(contData.strictMin)];
figure
title(['$v_\mathrm{avg}=$',num2str(vAVG)],'Interpreter','latex')
hold on
grid on
for i=1:(numel(idx)-1)
    range = idx(i):idx(i+1);
    if contData.strictMin(idx(i)+1)>0
        plot(contData.gamma(range)*180/pi,contData.T(range),'b')
    else
        plot(contData.gamma(range)*180/pi,contData.T(range),'b--')
    end
end
xlabel('$\gamma~[^\circ]$','Interpreter','latex')
ylabel('$T$','Interpreter','latex')

figure
title(['$v_\mathrm{avg}=$',num2str(vAVG)],'Interpreter','latex')
hold on
grid on
for i=1:(numel(idx)-1)
    range = idx(i):idx(i+1);
    if contData.strictMin(idx(i)+1)>0
        plot(contData.gamma(range)*180/pi,contData.cost(range),'b')
    else
        plot(contData.gamma(range)*180/pi,contData.cost(range),'b--')
    end
end
xlabel('$\gamma~[^\circ]$','Interpreter','latex')
ylabel('cost','Interpreter','latex')

figure
title(['$v_\mathrm{avg}=$',num2str(vAVG)],'Interpreter','latex')
hold on
grid on
plot(repmat(contData.gamma*180/pi,[nXi,1])',contData.xi')
xlabel('$\gamma~[^\circ]$','Interpreter','latex')
ylabel('$\xi_i$','Interpreter','latex')

figure
title(['$v_\mathrm{avg}=$',num2str(vAVG)],'Interpreter','latex')
hold on
grid on
plot(repmat(contData.gamma*180/pi,[n+2,1])',contData.lambda')
xlabel('$\gamma~[^\circ]$','Interpreter','latex')
ylabel('$\lambda_i$','Interpreter','latex')

