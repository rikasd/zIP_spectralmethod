function [M,C,G] = getDynamics_nonlinDIP(lumped_params,q,Dq)
%GETDYNAMICS_NONLINDIP Get dynamics parameters for the 
%                      double-inverted-pendulum (DIP) model
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
% q  = [ q_ankle;  q_hip] ankle angle and hip angle
% Dq = [Dq_ankle; Dq_hip] ankle angular velocity and hip angular velocity
%
% Outputs:
% M = inertial forces
% C = coreolis and centrifugal forces
% G = gravitational forces
%
% Rika Sugimoto Dimitrova
% 2024-02-10
% Reference: Shiozawa et al. 2021, Appendices 2-4

g = 9.81;

m1 = lumped_params.m1;
m2 = lumped_params.m2;
c1 = lumped_params.c1;
c2 = lumped_params.c2;
j1 = lumped_params.j1;
j2 = lumped_params.j2;
L1 = lumped_params.L1;

M = [(j1+m1*c1^2)+(j2+m2*c2^2)+m2*(L1^2+2*L1*c2*cos(q(2))), ...
     (j2+m2*c2^2)+m2*L1*c2*cos(q(2));...
     (j2+m2*c2^2)+m2*L1*c2*cos(q(2)),...
     (j2+m2*c2^2)];
C = m2*L1*c2*sin(q(2))*[2*Dq(2), Dq(2); Dq(1), 0];
G = -g*[m1*c1*sin(q(1))+m2*(L1*sin(q(1))+c2*sin(q(1)+q(2)));...
        m2*c2*sin(q(1)+q(2))];

end