function [A_cl,B_cl,C_cl,D_cl] = getModel_linDIP(lumped_params,controller_params)
%GETMODEL_LINDIP Get full-state-feedback closed-loop (cl) state-space model 
%               of linearized double-inverted pendulum (DIP) with 
%               linear quadratic regulator (LQR) controller 
%
% Inputs:
% lumped_params = struct containing lumped parameters for DIP
%   'm1' : mass of link 1 (ankle to hip)
%   'm2' : mass of link 2 (above hip)
%   'c1' : distance from ankle joint to center of mass of link 1
%   'c2' : distance from ankle joint to center of mass of link 2
%   'j1' : moment of inertia of link 1 about its center of mass
%   'j2' : moment of inertia of link 2 about its center of mass
%   'L1' : length of link 1
% controller_params = struct containing controller parameters for the LQR
%   'Q'  : 4x4 state penalty matrix
%   'R'  : 2x2 control-input penalty matrix
%
% Outputs:
% State-space model matrices A_cl, B_cl, C_cl, D_cl
%   A = dynamics matrix
%   B = input matrix
%   C = output matrix
%   D = direct feedthrough matrix
%
% Rika Sugimoto Dimitrova
% 2024-02-10
% Reference: Shiozawa et al. 2021, Appendices 2-4

m1 = lumped_params.m1;
m2 = lumped_params.m2;

q_eq = [0;0]; Dq_eq = [0;0];
[M,~,~] = getDynamics_nonlinDIP(lumped_params,q_eq,Dq_eq);
[J_CoM,DJ_CoM,J_G] = getJacobians_nonlinDIP(lumped_params,q_eq,Dq_eq);

% open-loop system state-space matrices
A_ol = [zeros(2) eye(2); -M\J_G zeros(2)];
B_ol = [zeros(2); M\eye(2)];

C_1 = [0 0 0 0];
D_1 = [1 0];

J_CoM_x = J_CoM(1,:);
DJ_CoM_x = DJ_CoM(1,:);
J_2 = -(m1+m2)*[DJ_CoM_x J_CoM_x];
C_2 = J_2*A_ol;
D_2 = J_2*B_ol;
C_ol = [C_1; C_2];
D_ol = [D_1; D_2];

% LQR controller gains
K_x = lqr(A_ol,B_ol,controller_params.Q,controller_params.R); 

A_cl = A_ol - B_ol*K_x;
B_cl = B_ol;
C_cl = C_ol - D_ol*K_x;
D_cl = D_ol;

end

