%% set up
addpath('utilities');
addpath('dynamics');
addpath('animations');

n   = 4; % state dimension

% passive solution at vAVG=0.1
load passiveGait_01.mat
vAVG      = passiveGait_01.vAVG;
T_pas     = passiveGait_01.T;
x0_pas    = passiveGait_01.x0;
gamma_pas = passiveGait_01.gamma;

gamma_end = 0;

fDstepSize = 1e-9; % finite difference step size
rootFunctionTolerance = 1e-8;
odeOpts = odeset('RelTol',1e-9,'AbsTol',1e-10); % ode solver
h = 1e-2; % fixed step-size in continuation

% reconstruct multipliers from passive solution
p0_pas     = zeros(4,1);
q_pas      = 1/(vAVG*T_pas);
lambda_pas = zeros(2,1);

u0_pas = 0;

% X = [T,x0,p0,q,u0,lambda,gamma,vAVG]
% X = [z;vAVG]
X0 = [T_pas;
      x0_pas;
      p0_pas;
      q_pas;
      u0_pas;
      lambda_pas;
      gamma_pas;
      vAVG];

rFun = @(z) resIndirect([z;vAVG],n,odeOpts); % fix vAVG
% z = [T,x0,p0,q,u0,lambda,gamma]
z_init = X0(1:end-1);

contData = indirectContinuation(z_init,rFun,gamma_end,n,...
                                fDstepSize,rootFunctionTolerance,h);

contData.gamma = contData.sigma;
contData.vAVG  = vAVG*ones(size(contData.T));
contData = rmfield(contData, 'sigma');

%% animations
% find passive gaits
[~,idxs] = sort(abs(contData.u0));
% passive gait with T_long
idx = idxs(1); % idx=1
[t_,x_]=getTrajectoriesIndirect(contData.T(idx),...
                                contData.x0(:,idx),...
                                contData.p0(:,idx),...
                                contData.q(idx),...
                                contData.u0(:,idx),...
                                contData.gamma(idx),...
                                contData.vAVG(idx),...
                                odeOpts);
getAnimationCompassGait(interp1(t_,x_,linspace(0,t_(end),50)'),0,1,'passive_T_long')

% passive gait with T_short
idx = idxs(2);
[t_,x_]=getTrajectoriesIndirect(contData.T(idx),...
                                contData.x0(:,idx),...
                                contData.p0(:,idx),...
                                contData.q(idx),...
                                contData.u0(:,idx),...
                                contData.gamma(idx),...
                                contData.vAVG(idx),...
                                odeOpts);
getAnimationCompassGait(interp1(t_,x_,linspace(0,t_(end),50)'),0,1,'passive_T_short')

% level ground walking
idx = idxs(end);
[t_,x_]=getTrajectoriesIndirect(contData.T(idx),...
                                contData.x0(:,idx),...
                                contData.p0(:,idx),...
                                contData.q(idx),...
                                contData.u0(:,idx),...
                                contData.gamma(idx),...
                                contData.vAVG(idx),...
                                odeOpts);
getAnimationCompassGait(interp1(t_,x_,linspace(0,t_(end),50)'),0,1,'level_ground')


%% plots
figure
title(['$v_\mathrm{avg}=$',num2str(vAVG)],'Interpreter','latex')
hold on
grid on
plot(contData.gamma*180/pi,contData.T)
xlabel('$\gamma~[^\circ]$','Interpreter','latex')
ylabel('$T$','Interpreter','latex')

figure
title(['$v_\mathrm{avg}=$',num2str(vAVG)],'Interpreter','latex')
hold on
grid on
plot(contData.gamma*180/pi,contData.cost)
xlabel('$\gamma~[^\circ]$','Interpreter','latex')
ylabel('cost','Interpreter','latex')

figure
title(['$v_\mathrm{avg}=$',num2str(vAVG)],'Interpreter','latex')
hold on
grid on
plot(repmat(contData.gamma*180/pi,[n,1])',contData.p0')
xlabel('$\gamma~[^\circ]$','Interpreter','latex')
ylabel('$p_0$','Interpreter','latex')

figure
title(['$v_\mathrm{avg}=$',num2str(vAVG)],'Interpreter','latex')
hold on
grid on
plot(repmat(contData.gamma*180/pi,[n,1])',contData.q')
xlabel('$\gamma~[^\circ]$','Interpreter','latex')
ylabel('$q$','Interpreter','latex')

figure
title(['$v_\mathrm{avg}=$',num2str(vAVG)],'Interpreter','latex')
hold on
grid on
plot(repmat(contData.gamma*180/pi,[2,1])',contData.lambda')
xlabel('$\gamma~[^\circ]$','Interpreter','latex')
ylabel('$\lambda_i$','Interpreter','latex')
