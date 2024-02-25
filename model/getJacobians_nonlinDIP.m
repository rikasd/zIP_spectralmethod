function [J_CoM,DJ_CoM,J_G] = getJacobians_nonlinDIP(lumped_params,q,Dq)
%GETJACOBIANS_NONLINDIP Get jacobians of center of mass and gravity force 
%                       of double-inverted-pendulum (DIP) model
%
% Reference: Shiozawa et al. 2021, Appendices 2-4
%
% Inputs:
% lumped_params = struct containing lumped parameters for DIP
%   'm1' : mass of link 1 (ankle to hip)
%   'm2' : mass of link 2 (above hip)
%   'c1' : distance from ankle joint to center of mass of link 1
%   'c2' : distance from ankle joint to center of mass of link 2
%   'L1' : length of link 1
%
% Outputs:
% J_CoM  = Jacobian of the center of mass (CoM) position with respect to
%          the ankle and jip joint angles q = [q_ankle; q_hip]
% DJ_CoM = Time derivative of J_CoM
% J_G    = Jacobian of gravitational forces G
%
% Rika Sugimoto Dimitrova
% 2024-02-10

g = 9.81;

m1 = lumped_params.m1;
m2 = lumped_params.m2;
c1 = lumped_params.c1;
c2 = lumped_params.c2;
L1 = lumped_params.L1;

M1 = m1/(m1+m2);
M2 = m2/(m1+m2);

J_CoM1 = [M1*c1*cos(q(1))+M2*(L1*cos(q(1))+c2*cos(q(1)+q(2)));...
          M1*c1*sin(q(1))+M2*(L1*sin(q(1))+c2*sin(q(1)+q(2)))];
J_CoM2 = [M2*c2*cos(q(1)+q(2)); M2*c2*sin(q(1)+q(2))];
J_CoM = -[J_CoM1 J_CoM2];

DJ_CoM = [M1*c1*sin(q(1))*Dq(1)+M2*(L1*sin(q(1))*Dq(1)+c2*sin(q(1)+q(2))*(Dq(1)+Dq(2))) ...
          M2*c2*sin(q(1)+q(2))*(Dq(1)+Dq(2));...
         -M1*c1*cos(q(1))*Dq(1)-M2*(L1*cos(q(1))*Dq(1)+c2*cos(q(1)+q(2))*(Dq(1)+Dq(2))) ...
         -M2*c2*cos(q(1)+q(2))*(Dq(1)+Dq(2))];

J_G = -g*[m1*c1*cos(q(1))+m2*(L1*cos(q(1))+c2*cos(q(1)+q(2))), m2*c2*cos(q(1)+q(2));
          m2*c2*cos(q(1)+q(2)), m2*c2*cos(q(1)+q(2))];

end