%% Symbolic Derivation of Equations of Motion (EOM)
% This script derives the Equations of Motion (EOM) using the Euler-Lagrange
% method for a compass-gait walker, represented in floating base coordinates. 
% It also includes the transformation of the EOM into minimal coordinates and 
% the computation of their derivatives. The dynamics are automatically generated 
% and exported into MATLAB functions.

%% Floating Base Variables
% Number of floating base variables and their symbolic definition
n = 4;  % Number of floating base variables
syms x y sw st real   % Generalized coordinates
syms dx dy dsw dst real % Generalized velocities
syms u real  % Input variable (control)

q  = [x y sw st];  % Generalized coordinates vector
dq = [dx dy dsw dst];  % Generalized velocities vector

gamma = sym('gamma', 'real');  % Slope of the ground (incline)

%% Model Parameters
% Define symbolic variables for the model's physical parameters
m   = sym('m', 'real');    % Mass of the leg
m_h = sym('m_h', 'real');  % Mass of the hip
g   = sym('g', 'real');    % Gravitational acceleration
a   = sym('a', 'real');    % Length parameter a (from hip to center of mass)
b   = sym('b', 'real');    % Length parameter b (from hip to foot)

l = a + b;  % Total length of the leg

% Gravitational force vector in terms of incline gamma
g_vec = [sin(gamma); -cos(gamma)] * g;

%% --- Kinematics ---
% Centers of Gravity (CoG) and velocities computed via Jacobians

% Positions of the centers of gravity (CoG)
CoG_T   = [x; y];  % CoG of the torso (hip)
pos_Hip = CoG_T;  % Position of the hip (same as CoG of torso)
CoG_sw  = pos_Hip + b * [sin(sw); -cos(sw)];  % CoG of the swinging leg
CoG_st  = pos_Hip + b * [sin(st); -cos(st)];  % CoG of the stance leg

% Velocities of the centers of gravity
d_CoG_T  = jacobian(CoG_T, q) * dq.';  % Velocity of the torso
d_CoG_sw = jacobian(CoG_sw, q) * dq.';  % Velocity of the swinging leg
d_CoG_st = jacobian(CoG_st, q) * dq.';  % Velocity of the stance leg

%% --- Energies ---
% Define potential and kinetic energy expressions

% Potential Energy (due to gravity)
V = -m_h * g_vec' * CoG_T - m * g_vec' * CoG_sw - m * g_vec' * CoG_st;

% Kinetic Energy
T = 0.5 * ( m_h * sum(d_CoG_T.^2) + ...
            m * sum(d_CoG_sw.^2) + ...
            m * sum(d_CoG_st.^2) );

%% --- Euler-Lagrange Equations ---
% M * ddq + C * dq + G = 0
% Derive mass matrix (M), Coriolis matrix (C), and gravity vector (G) using
% Euler-Lagrange formalism.
[M, C, CMat, G] = eulerLagrange(T, V, q, dq); 

%% --- Contact Projection Matrices ---
% Compute the contact matrices for the stance and swing legs, and their derivatives

% Contact points of stance and swing legs
g_st = [x + l * sin(st); y - l * cos(st)];  % Stance leg
g_sw = [x + l * sin(sw); y - l * cos(sw)];  % Swing leg

% Compute contact matrices (W) and their time derivatives (W_dot)
[W_st, W_st_dot] = computeContactMatrix(g_st, q, dq);  % Stance leg contact matrix
[W_sw, W_sw_dot] = computeContactMatrix(g_sw, q, dq);  % Swing leg contact matrix

% Compute discrete map for collision dynamics
Gd    = W_sw' * (M \ W_sw); 
Delta = blkdiag(eye(2), [0,1;1,0]) * (eye(4) - M \ (W_sw * (Gd \ W_sw')));  % Discrete map for impact

%% Transformation into Minimal Coordinates
% Transform the dynamics into minimal coordinates z = [swM, stM]

syms swM stM real  % Minimal coordinates
syms dswM dstM real  % Derivatives of minimal coordinates
z  = [swM; stM];    % Minimal coordinates vector
dz = [dswM; dstM];  % Velocities in minimal coordinates

% Define transformation from floating base coordinates to minimal coordinates
q_z = [-l * sin(stM + gamma); l * cos(stM + gamma); swM + gamma; stM + gamma];
BTrafo = jacobian(q_z, z);  % Transformation matrix

% Transform velocities into minimal coordinates
dqNew = BTrafo * dz;

% Compute time derivative of transformation matrix
BdtTrafo = BTrafo;
for i = 1:size(BdtTrafo,1)
    for j = 1:size(BdtTrafo,2)
        BdtTrafo(i,j) = jacobian(BdtTrafo(i,j), z) * dz;
    end
end

% Minimal coordinate mass matrix and Coriolis term
M_min = simplify(BTrafo' * subs(M, q', q_z) * BTrafo);
c_min = simplify(BTrafo' * subs(C, [q'; dq'], [q_z; dqNew]) + BTrafo' * subs(M, q', q_z) * BdtTrafo * dz);

% Potential energy and control input in minimal coordinates
G_min = simplify(BTrafo' * subs(G, q', q_z));
B_min = simplify(BTrafo' * [0 0 -1 1]');

% Minimal coordinate system dynamics
f_min = [dz; simplify(M_min \ (B_min * u - G_min - c_min))];

% Minimal coordinate discrete map
Delta_min = simplify(expand([zeros(2), eye(2)] * subs(Delta, q', q_z) * BTrafo));
g_min = [stM; swM; Delta_min * dz];

%% Automatic Function Generation
% Generate MATLAB functions for the dynamics and their derivatives

% Substitute numerical values for parameters in f_min and g_min
f = subs(f_min, {a, b, g, m, m_h}, {0.5, 0.5, 1, 0.25, 0.5});
g = subs(g_min, {a, b, g, m, m_h}, {0.5, 0.5, 1, 0.25, 0.5});

x = [z; dz];  % State vector
dgdx = jacobian(g, x);  % Jacobian of g with respect to x
dfdx = jacobian(f, x);  % Jacobian of f with respect to x
dfdu = jacobian(f, u);  % Jacobian of f with respect to u

% Generate MATLAB functions and export them to files
matlabFunction(f, 'File', 'fAUTO', 'Vars', {x, u});
matlabFunction(dfdx, 'File', 'f_xAUTO', 'Vars', {x, u});
matlabFunction(dfdu, 'File', 'f_uAUTO', 'Vars', {x, u});
matlabFunction(g, 'File', 'gAUTO', 'Vars', {x});
matlabFunction(dgdx, 'File', 'g_xAUTO', 'Vars', {x});

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