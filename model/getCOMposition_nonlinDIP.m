function [p_O_CoM] = getCOMposition_nonlinDIP(lumped_params,q)
%GETCOMPOSITION_NONLINDIP Compute current position of center of mass (CoM)
%
% Inputs:
% lumped_params = struct containing lumped parameters for DIP
%   'm1' : mass of link 1 (ankle to hip)
%   'm2' : mass of link 2 (above hip)
%   'c1' : distance from ankle joint to center of mass of link 1
%   'c2' : distance from ankle joint to center of mass of link 2
%   'L1' : length of link 1
% q  = [ q_ankle;  q_hip] ankle angle and hip angle
%
% Outputs:
% p_O_CoM = [x_CoM; z_CoM] CoM horizontal and vertical coordinates 
%
% Rika Sugimoto Dimitrova
% 2024-02-11

m1 = lumped_params.m1;
m2 = lumped_params.m2;
c1 = lumped_params.c1;
c2 = lumped_params.c2;
L1 = lumped_params.L1;

M1 = m1/(m1+m2);
M2 = m2/(m1+m2);

x_CoM = -M1*c1*sin(q(1)) - M2*(L1*sin(q(1)) + c2*sin(q(1)+q(2)));
z_CoM =  M1*c1*cos(q(1)) + M2*(L1*cos(q(1)) + c2*cos(q(1)+q(2)));

p_O_CoM = [x_CoM; z_CoM];

end
