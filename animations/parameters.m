function p = parameters()
% PARAMETERS Computes the model parameters for the simplest walking model.
%
% P = PARAMETERS() returns a struct P with fields:
%       - P.m: the wheel's total mass
%       - P.L: the wheel's leg length
%       - P.g: gravity
%       - P.r_g: radius of gyration
%       - P.alpha: leg angle
%       - P.gamma: slope of the ramp

%% model parameters
p.L = 1; % leg length [units don't matter, could be m]
end