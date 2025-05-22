%% Symbolic Derivation of Equations of Motion (EOM)
% This script derives the Equations of Motion (EOM) using the Euler-Lagrange
% method for the robot RABBIT, represented in floating base coordinates. 
% It also includes the transformation of the EOM into minimal coordinates and 
% the computation of their derivatives. The dynamics are automatically generated 
% and exported into MATLAB functions.

%% Floating Base Variables
% Number of floating base variables and their symbolic definition
n = 7;  % Number of floating base variables
syms x y q1 q2 q3 q4 q5 real  % Generalized coordinates
syms dx dy dq1 dq2 dq3 dq4 dq5 real % Generalized velocities
syms u1 u2 u3 u4 real  % Input variable (control)

q  = [q1; q2; q3; q4; q5; x; y];         % Generalized coordinates vector
dq = [dq1; dq2; dq3; dq4; dq5; dx; dy];  % Generalized velocities vector
u  = [u1; u2; u3; u4];                   % Generalized input vector

gamma = sym('gamma', 'real');  % Slope of the ground (incline)

%% Model Parameters
% Define symbolic variables for the model's physical parameters
M_T     = sym('M_T',{'real','positive'}); %_T: torso
M_f     = sym('M_f',{'real','positive'}); %_f: femur
M_t     = sym('M_t',{'real','positive'}); %_t: tibia
l_T     = sym('l_T',{'real','positive'});
l_f     = sym('l_f',{'real','positive'});
l_t     = sym('l_t',{'real','positive'});
I_T     = sym('I_T',{'real','positive'});
I_f     = sym('I_f',{'real','positive'});
I_t     = sym('I_t',{'real','positive'});
p_T     = sym('p_T',{'real','positive'});
p_f     = sym('p_f',{'real','positive'});
p_t     = sym('p_t',{'real','positive'});
g       = sym('g',{'real','positive'});

% Gravitational force vector in terms of incline gamma
g_vec = [sin(gamma); -cos(gamma)] * g;

%% --- Kinematics ---
% Centers of Gravity (CoG) and velocities computed via Jacobians

% Positions of the centers of gravity (CoG)
CoG_T       = [x; y];  % CoG of the torso 
Angle_T     = q5;
pos_Hip     = CoG_T+p_T*[sin(Angle_T);-cos(Angle_T)];  % Position of the hip
Angle_f_sw  = q2+q5;
Angle_f_st  = q1+q5;
CoG_f_sw    = pos_Hip + p_f * [-sin(Angle_f_sw); cos(Angle_f_sw)];  % CoG of the swinging femur
CoG_f_st    = pos_Hip + p_f * [-sin(Angle_f_st); cos(Angle_f_st)];  % CoG of the stance femur
pos_Knee_sw = pos_Hip + l_f * [-sin(Angle_f_sw); cos(Angle_f_sw)];  % Position of the swing knee
pos_Knee_st = pos_Hip + l_f * [-sin(Angle_f_st); cos(Angle_f_st)];  % Position of the stance knee
Angle_t_sw  = q2+q4+q5;
Angle_t_st  = q1+q3+q5;
CoG_t_sw    = pos_Knee_sw + p_t * [-sin(Angle_t_sw); cos(Angle_t_sw)];  % CoG of the swinging tibia
CoG_t_st    = pos_Knee_st + p_t * [-sin(Angle_t_st); cos(Angle_t_st)];  % CoG of the stance tibia
pos_Foot_sw = pos_Knee_sw + l_t * [-sin(Angle_t_sw); cos(Angle_t_sw)];  % Position of the swing foot
pos_Foot_st = pos_Knee_st + l_t * [-sin(Angle_t_st); cos(Angle_t_st)];  % Position of the stance foot

% Linear and Angular Velocities of the centers of gravity
d_CoG_T  = jacobian(CoG_T, q) * dq;   % Velocity of the torso
d_CoG_f_sw = jacobian(CoG_f_sw, q) * dq;  % Velocity of the swinging femur
d_CoG_f_st = jacobian(CoG_f_st, q) * dq;  % Velocity of the stance femur
d_CoG_t_sw = jacobian(CoG_t_sw, q) * dq;  % Velocity of the swinging tibia
d_CoG_t_st = jacobian(CoG_t_st, q) * dq;  % Velocity of the stance tibia

d_Angle_T = jacobian(Angle_T, q) * dq;
d_Angle_f_sw = jacobian(Angle_f_sw, q) * dq;
d_Angle_f_st = jacobian(Angle_f_st, q) * dq;
d_Angle_t_sw = jacobian(Angle_t_sw, q) * dq;
d_Angle_t_st = jacobian(Angle_t_st, q) * dq;

%% --- Energies ---
% Define potential and kinetic energy expressions

% Potential Energy (due to gravity)
V = -M_T*g_vec'*CoG_T ...
    -M_f*g_vec'*CoG_f_sw ...
    -M_f*g_vec'*CoG_f_st ...
    -M_t*g_vec'*CoG_t_sw ...
    -M_t*g_vec'*CoG_t_st;

% Kinetic Energy
T = 0.5 * ( + M_T * sum(d_CoG_T.^2) ...
            + M_f * sum(d_CoG_f_sw.^2) ...
            + M_f * sum(d_CoG_f_st.^2) ...
            + M_t * sum(d_CoG_t_sw.^2) ...
            + M_t * sum(d_CoG_t_st.^2) ...
            + I_T * sum(d_Angle_T.^2) ...
            + I_f * sum(d_Angle_f_sw.^2) ...
            + I_f * sum(d_Angle_f_st.^2) ...
            + I_t * sum(d_Angle_t_sw.^2) ...
            + I_t * sum(d_Angle_t_st.^2) );

%% --- Euler-Lagrange Equations ---
% M * ddq + C * dq + G = 0
% Derive mass matrix (M), Coriolis matrix (C), and gravity vector (G) using
% Euler-Lagrange formalism.
[M, C, CMat, G] = eulerLagrange(T, V, q', dq'); 

%% --- Contact Projection Matrices ---
% Compute the contact matrices for the stance and swing legs, and their derivatives

% - Contact Forces -
g_st = pos_Foot_st;                           
g_sw = pos_Foot_sw;                   

% Compute contact matrices (W) and their time derivatives (W_dot)
[W_st, W_st_dot] = computeContactMatrix(g_st, q', dq');  % Stance leg contact matrix
[W_sw, W_sw_dot] = computeContactMatrix(g_sw, q', dq');  % Swing leg contact matrix

%% Transformation into Minimal Coordinates
% Transform the dynamics into minimal coordinates z = [q1, q2, q3, q4, q5]
syms q1M  q2M  q3M  q4M  q5M  real  % Minimal coordinates
syms dq1M dq2M dq3M dq4M dq5M real  % Derivatives of minimal coordinates
z  = [q1M;  q2M;  q3M;  q4M;  q5M];   % Minimal coordinates vector
dz = [dq1M; dq2M; dq3M; dq4M; dq5M];  % Velocities in minimal coordinates

% Define transformation from floating base coordinates to minimal coordinates
q_z = [ q1M; ...
        q2M; ...
        q3M; ...
        q4M; ...
        q5M; ...
        + l_f*sin(q1M+q5M) - p_T*sin(q5M) + l_t*sin(q1M+q3M+q5M);...
        - l_f*cos(q1M+q5M) + p_T*cos(q5M) - l_t*cos(q1M+q3M+q5M); ...
       ];
BTrafo = jacobian(q_z, z);  % Transformation matrix

% Transform velocities into minimal coordinates
dq_z = BTrafo * dz;

% M(q) -> M(z)
M = subs(M, q, q_z);

% W_sw(q) -> W_sw(z)
W_sw = subs(W_sw, q, q_z);

% Compute time derivative of transformation matrix
BdtTrafo = BTrafo;
for i = 1:size(BdtTrafo,1)
    for j = 1:size(BdtTrafo,2)
        BdtTrafo(i,j) = jacobian(BdtTrafo(i,j), z) * dz;
    end
end

% Minimal coordinate mass matrix and Coriolis term
M_min = simplify(BTrafo' * M * BTrafo);
c_min = simplify(BTrafo' * subs(C, [q; dq], [q_z; dq_z]) + BTrafo' * M * BdtTrafo * dz);

% Potential energy and control input in minimal coordinates
G_min = simplify(BTrafo' * subs(G, q, q_z));
B_min = simplify(BTrafo' * [eye(4);zeros(3,4)]);

% Minimal coordinate system dynamics
%f_min = [dz; simplify(M_min \ (B_min * u - G_min - c_min))];

% Minimal coordinate discrete map
%Gd    = W_sw' * (M \ W_sw); 
%Delta = blkdiag([0 1 0 0; 1 0 0 0; 0 0 0 1; 0 0 1 0],eye(3)) * (eye(7) - M \ (W_sw * (Gd \ W_sw')));  % Discrete map for impact
%Delta_min = simplify(expand([eye(5), zeros(5,2)] * Delta * BTrafo));
%g_min = [q2M; q1M; q4M; q3M; q5M; Delta_min * dz];


%% Automatic Function Generation
% Generate MATLAB functions for the dynamics and their derivatives

% Substitute numerical values for parameters (normalized)
m0_val  = 12+2*(6.8+3.2); % total Mass
l0_val = 0.8; % leg length
g0_val  = 9.81; % gravity
I0_val  = m0_val*l0_val^2; % inertia

M_T_val = 12/m0_val;
M_f_val = 6.8/m0_val;
M_t_val = 3.2/m0_val;

l_T_val = 0.63/l0_val;
l_f_val = 0.4/l0_val;
l_t_val = 0.4/l0_val;

I_T_val = 1.33/I0_val;
I_f_val = 0.47/I0_val;
I_t_val = 0.2/I0_val;

p_T_val = 0.24/l0_val;
p_f_val = 0.11/l0_val;
p_t_val = 0.24/l0_val;

g_val = 1; % 9.81/g0_val

M = subs(M, {g, M_T,M_f,M_t, l_T,l_f,l_t, I_T,I_f,I_t, p_T,p_f,p_t}, {g_val, M_T_val,M_f_val,M_t_val, l_T_val,l_f_val,l_t_val, I_T_val,I_f_val,I_t_val, p_T_val,p_f_val,p_t_val});
W_sw = subs(W_sw, {g, M_T,M_f,M_t, l_T,l_f,l_t, I_T,I_f,I_t, p_T,p_f,p_t}, {g_val, M_T_val,M_f_val,M_t_val, l_T_val,l_f_val,l_t_val, I_T_val,I_f_val,I_t_val, p_T_val,p_f_val,p_t_val});
BTrafo = subs(BTrafo, {g, M_T,M_f,M_t, l_T,l_f,l_t, I_T,I_f,I_t, p_T,p_f,p_t}, {g_val, M_T_val,M_f_val,M_t_val, l_T_val,l_f_val,l_t_val, I_T_val,I_f_val,I_t_val, p_T_val,p_f_val,p_t_val});

M_min = subs(M_min, {g, M_T,M_f,M_t, l_T,l_f,l_t, I_T,I_f,I_t, p_T,p_f,p_t}, {g_val, M_T_val,M_f_val,M_t_val, l_T_val,l_f_val,l_t_val, I_T_val,I_f_val,I_t_val, p_T_val,p_f_val,p_t_val});
c_min = subs(c_min, {g, M_T,M_f,M_t, l_T,l_f,l_t, I_T,I_f,I_t, p_T,p_f,p_t}, {g_val, M_T_val,M_f_val,M_t_val, l_T_val,l_f_val,l_t_val, I_T_val,I_f_val,I_t_val, p_T_val,p_f_val,p_t_val});
G_min = subs(G_min, {g, M_T,M_f,M_t, l_T,l_f,l_t, I_T,I_f,I_t, p_T,p_f,p_t}, {g_val, M_T_val,M_f_val,M_t_val, l_T_val,l_f_val,l_t_val, I_T_val,I_f_val,I_t_val, p_T_val,p_f_val,p_t_val});

x = [z; dz];  % State vector
dMdx = sym(zeros(7,7,10)); % Jacobian of M with respect to x
for i = 1:10
    dMdx(:,:,i) = diff(M, x(i));
end
dW_swdx = sym(zeros(7,2,10)); % Jacobian of W_sw with respect to x
for i = 1:10
    dW_swdx(:,:,i) = diff(W_sw, x(i));
end
dBTrafodx = sym(zeros(7,5,10)); % Jacobian of BTrafo with respect to x
for i = 1:10
    dBTrafodx(:,:,i) = diff(BTrafo, x(i));
end
dM_mindx = sym(zeros(5,5,10)); % Jacobian of M_min with respect to x
for i = 1:10
    dM_mindx(:,:,i) = diff(M_min, x(i));
end
dc_mindx = jacobian(c_min, x); % Jacobian of c_min with respect to x
dG_mindx = jacobian(G_min, x); % Jacobian of c_min with respect to x
dc_mindgamma = diff(c_min, gamma); % Jacobian of c_min with respect to gamma
dG_mindgamma = diff(G_min, gamma); % Jacobian of c_min with respect to gamma


% Generate MATLAB functions and export them to files
matlabFunction(M, 'File', 'MAUTO', 'Vars', {x});
matlabFunction(W_sw, 'File', 'W_swAUTO', 'Vars', {x});
matlabFunction(BTrafo, 'File', 'BTrafoAUTO', 'Vars', {x});
matlabFunction(M_min, 'File', 'M_minAUTO', 'Vars', {x});
matlabFunction(c_min, 'File', 'c_minAUTO', 'Vars', {x,gamma});
matlabFunction(G_min, 'File', 'G_minAUTO', 'Vars', {x,gamma});

matlabFunction(dMdx, 'File', 'dMdxAUTO', 'Vars', {x});
matlabFunction(dW_swdx, 'File', 'dW_swdxAUTO', 'Vars', {x});
matlabFunction(dBTrafodx, 'File', 'dBTrafodxAUTO', 'Vars', {x});
matlabFunction(dM_mindx, 'File', 'dM_mindxAUTO', 'Vars', {x});
matlabFunction(dc_mindx, 'File', 'dc_mindxAUTO', 'Vars', {x,gamma});
matlabFunction(dG_mindx, 'File', 'dG_mindxAUTO', 'Vars', {x,gamma});
matlabFunction(dc_mindgamma, 'File', 'dc_mindgammaAUTO', 'Vars', {x,gamma});
matlabFunction(dG_mindgamma, 'File', 'dG_mindgammaAUTO', 'Vars', {x,gamma});


%% Helper Functions

function [M, C, CMat, G] = eulerLagrange(T, V, q, dq)
    % Compute the mass matrix (M), Coriolis forces (C), and gravity vector (G)
    % using the Euler-Lagrange formalism.
    
    % Lagrangian (T - V)
    L = simplify(T - V);

    % Partial derivatives
    dT_dq   = jacobian(T, q).';  % Partial derivative of kinetic energy w.r.t. q
    dV_dq   = jacobian(V, q).';  % Partial derivative of potential energy w.r.t. q
    dL_dqdt = jacobian(L, dq).'; % Partial derivative of Lagrangian w.r.t. dq

    % Time derivative of dL/dqdt
    dd_L_dqdt2  = jacobian(dL_dqdt, dq);
    d_dLdqdt_dq = jacobian(dL_dqdt, q);

    % Assign matrices
    M = simplify(dd_L_dqdt2);  % Mass matrix
    CMat = -simplify(expand(jacobian(0.5 * dq * M, q).' - d_dLdqdt_dq));  % Coriolis matrix (expanded form)
    C = -simplify(expand(dT_dq - d_dLdqdt_dq * dq.'));  % Coriolis vector
    G = simplify(expand(dV_dq));  % Gravity vector
end

function [W, W_dot] = computeContactMatrix(con, q, dq)
    % Compute the contact projection matrix W and its time derivative W_dot
    % for a given constraint 'con'.
    
    W = jacobian(con, q)';  % Contact matrix
    W_dot = W;

    % Compute time derivative of each column of W
    for j = 1:size(W, 2)
        W_dot(:, j) = jacobian(W(:, j), q) * dq';  % Time derivative of contact matrix
    end
end