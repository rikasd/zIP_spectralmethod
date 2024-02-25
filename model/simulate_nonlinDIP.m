function [output_struct] = simulate_nonlinDIP(input_struct)
%SIMULATE_NONLINDIP Simulate nonlinear double-inverted pendulum (DIP) model
%
% Inputs:
% input_struct = struct of simulation parameters, including:
%   'lumped_params'     : struct of lumped parameters for DIP
%     'm1' : mass of link 1 (ankle to hip)
%     'm2' : mass of link 2 (above hip)
%     'c1' : distance from ankle joint to center of mass of link 1
%     'c2' : distance from ankle joint to center of mass of link 2
%     'j1' : moment of inertia of link 1 about its center of mass
%     'j2' : moment of inertia of link 2 about its center of mass
%     'L1' : length of link 1
%   'simFreq_Hz'        : simulation frequency (Hz), default 1000 Hz
%   'sampFreq_Hz'       : output data sampling frequency (Hz),
%                         default 100 Hz
%   'simDuration_s'     : simulation duration (s), default 60 s
%   'motorNoiseLvL_Nm'  : standard deviation of actuator noise (Nm),
%                         default 10 Nm
%   'motorNoiseRatio'   : ankle actuator noise / hip actuator noise,
%                         default 1.6
%   'controller_params' : struct of LQR cost function weighting matrices
%     'Q'  : 4x4 state penalty matrix, default ones(4)
%     'R'  : 2x2 control-input penalty matrix, default 1e6*[0.3 0; 0 1/0.3]
%
% Outputs:
% output_struct = struct of simulation output, including:
%   'COP'         : horizontal coordinate of center of pressure (COP)
%   'COM'         : [COM_x; COM_z] horizontal and vertical coordinates of
%                                  center of mass (COM)
%   'FootForce'   : [F_x; F_z] horizontal and vertical components of
%                              the foot-force vector
%   'sampFreq_Hz' : data sampling frequency
%   'thetaF'      : foot-force vector orientation
%
% Rika Sugimoto-Dimitrova (rikasd@mit.edu)
% 2024-02-10
% Adapted code from Jongwoo Lee, PhD. and Kaymie Shiozawa, MS. (Shiozawa et al. 2021)

%% Model Parameters
if nargin ~=0
    lumped_params     = input_struct.lumped_params;
    simFreq_Hz        = input_struct.simFreq_Hz;
    sampFreq_Hz       = input_struct.sampFreq_Hz;
    t_0               = 0; 
    t_f               = input_struct.simDuration_s;
    motorNoiseLvL     = input_struct.motorNoiseLvL_Nm;
    motorNoiseRatio   = input_struct.motorNoiseRatio;
    noise_type        = input_struct.noise_type;
    controller_params = input_struct.controller_params;
else
    % Default parameters:
    %-- DIP model parameters --%
    lumped_params = ... % Lumped model parameters from Shiozawa et al. 2021
        struct('m1',26.30,'m2',42.88,...
               'c1',0.589,'c2',0.332,...
               'j1',1.400,'j2',2.227,...
               'L1',0.867,'L_COM',0.967);
    %-- Simulation setup --%
    simFreq_Hz      = 1000;
    sampFreq_Hz     = 100; 
    t_0             = 0;
    t_f             = 60;
    %-- Noise parameters --%
    motorNoiseLvL   = 0.01;
    motorNoiseRatio = 1.6;
    noise_type      = 'w';
    %-- LQR cost function weighting matrices --%
    gamma           = 1; 
    alpha           = 1e6;
    beta            = 0.3;
    controller_params.Q = gamma*eye(4); 
    controller_params.R = alpha*diag([beta, 1/beta]);
end

%% Simulation Setup
% initial condition
% q: joint ('relative, not spatial') angles
% q=[q_ankle, q_hip]
q_0 = [0;0];
Dq_0 = [0;0];

n_q = length(q_0); % number of joints

% x: state = [q;Dq]
x_0 = [q_0; Dq_0]; % initial state
n_x = length(x_0); % number of states

%% Run simulation
dt = 1/simFreq_Hz;    % timestep
t = t_0:dt:t_f; % time
N = length(t);  % sample number

% memory allocation 
x = zeros(n_x, N);         % 4xN matrix - trajectory 
Dx = zeros(n_x, N);
torque = zeros(n_q, N);
y = cell(1, N);            % 1x Nt structures - output 

% Uncorrelated Gaussian noise
motorNoise_white = motorNoiseLvL*randn([n_q,N]);
motorNoise = motorNoise_white;

if noise_type == 'l' % low-pass-filtered white noise (Peterka 2000 model)
    H = tf([1],[1 1/80]); % 80 s time constant
    Hd = c2d(H,dt);
    [num,den] = tfdata(Hd);
    motorNoise = filter(num{1},den{1},motorNoise_white')';
end

% run simulation
x_i = x_0; % initialize state

% Generate LQR controller gain
[~,~,K] = LQR_Controller(t_0, x_0, controller_params);

for i = 1:N
    t_i = t(i);
    u_i = LQR_Controller(t_i, x_i, controller_params);
    torque(:, i) = u_i;
    [Dx_i, y_i] = ...
        func_Dynamics(t_i, x_i, u_i, motorNoise(:,i));
    x(:, i) = x_i;
    Dx(:, i) = Dx_i;
    y{i} = y_i;
    
    % integrating
    % unwrapping states
    q_t = x_i(1:n_q, 1);
    Dq_t = x_i(n_q+1:end, 1);
    DDq_t = Dx_i(n_q+1:end, 1);
    % semi-implicit Euler integration
    Dq_next = Dq_t + dt*DDq_t;
    q_post = q_t + dt*Dq_next;
    x_i = [q_post; Dq_next];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% post process
COP_x = zeros(1,N); % horizontal CoP
COM_z = zeros(1,N); % vertical CoM
COM_x = zeros(1,N); % horizontal CoM
F_x   = zeros(1,N); % horizontal ground reactin force
F_z   = zeros(1,N); % vertical ground reaction force
% renaming component of output y from simulation above
count = 1;
for i = 1:N
    COP_x(:,count) = y{i}.COP;
    COM_z(:,count) = y{i}.COM_z;
    COM_x(:,count) = y{i}.COM_x; % 2021-06-21
    F_x(:,count)   = y{i}.Fx;
    F_z(:,count)   = y{i}.Fz;
    count = count + 1;
end

% Downsample data, if measurement frequency < simulaiton frequency
r = simFreq_Hz/sampFreq_Hz;
output_struct.COP = decimate(COP_x,r);
output_struct.COM = [decimate(COM_x,r); decimate(COM_z,r)];
output_struct.FootForce = [decimate(F_x,r); decimate(F_z,r)];
output_struct.sampFreq_Hz = sampFreq_Hz;
output_struct.thetaF = -decimate(F_x,r)./decimate(F_z,r);

%% sub function definitions
% Controller
% LQR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [u_LQR, A_cl, K_lqr] = LQR_Controller(t, x, param)
    Q = param.Q; R = param.R;
    if isempty(K) % finding K_lqr for u = -Kx
        [A_lin,B_lin] = getLinearDynamics();
        K_lqr = lqr(A_lin,B_lin,Q,R);
    else
        K_lqr = K;
    end
    [A_lin,B_lin] = getLinearDynamics();

    A_cl = A_lin-B_lin*K_lqr;
    
    q_eq = [0; 0]*pi/180;
    Dq_eq = [0; 0];
    
    x_eq = [q_eq; Dq_eq];
    u_LQR = -K_lqr*(x-x_eq);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plant Dynamics
function [Dx, output_struct] = ...
        func_Dynamics(t, x, u_ctrl, motorNoise)
    
    q = x(1:2,:);
    Dq = x(3:4,:);

    if isempty(motorNoise)
        motorNoise = 10*randn(length(q),1);
    end
    if isempty(motorNoiseRatio)
        motorNoiseRatio = 1.6; % ankle:hip noise ratio
    end
            
    % eom Mq'' + Cq' + G = u_ctrl + u_noise
    [M, C, G] = getDynamics_nonlinDIP(lumped_params,q,Dq);
    
    % input noise
    u_noise = motorNoise;
    u_noise(1) = motorNoiseRatio*u_noise(1);
    u = u_ctrl + u_noise;
    
    % solve accl 
    DDq = M\(-C*Dq -G + u);
    Dx = [Dq; DDq];
    
    % outputs
    if nargout > 1           
        % foot-ground interaction force
        F_vec = getFootForce_nonlinDIP(lumped_params,q,Dq,DDq);

        output_struct.Fx = F_vec(1);
        output_struct.Fz = F_vec(2);
                    
        % COP position
        output_struct.COP = u(1)/F_vec(2);
        
        % COM position
        p_O_CoM = getCOMposition_nonlinDIP(lumped_params,q);
        output_struct.COM_x = p_O_CoM(1);
        output_struct.COM_z = p_O_CoM(2);       
    end 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linearized Dynamics Matrices
function [A_lin,B_lin,J_G_eq] = getLinearDynamics()
    % equilibrium position
    q_eq = [0; 0]*pi/180; Dq_eq = [0;0];

    [M_eq, ~, ~] = getDynamics_nonlinDIP(lumped_params,q_eq,Dq_eq);
    [~, ~, J_G_eq]  = getJacobians_nonlinDIP(lumped_params,q_eq,Dq_eq);
    A_lin = [zeros(2),  eye(2); 
        -inv(M_eq)*J_G_eq, zeros(2)];
    B_lin = [zeros(2,2); inv(M_eq)];
end

end