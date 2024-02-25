% Double-inverted-pendulum lumped model parameters
lumped_params = ... % Lumped model parameters from Shiozawa et al. 2021
        struct('m1',26.30,'m2',42.88,...
               'c1',0.589,'c2',0.332,...
               'j1',1.400,'j2',2.227,...
               'L1',0.867,'L_COM',0.967);

% LQR cost function parameters
gamma = 1; alpha = 1e6; beta = 0.3; % from Shiozawa et al. 2021
controller_params.Q = gamma*eye(4); 
controller_params.R = alpha*diag([beta, 1/beta]);
    
% Simulation parameters
input_struct.lumped_params      = lumped_params;
input_struct.simFreq_Hz         = 1000;
input_struct.sampFreq_Hz        = 100;
input_struct.simDuration_s      = 60;
input_struct.motorNoiseLvL_Nm   = 0.1;
input_struct.motorNoiseRatio    = 0.9;
input_struct.noise_type         = 'w';
input_struct.controller_params  = controller_params;

N_trial = 40; % Number of simulations to run

zIP_params.method = 'cpsd';
zIP_params.window_size = 2^9;
zIP_params.detrend_option = 3; % linear detrend
zIP_params.window_taper_option = 2; % Hann
zIP_params.f_int = 0.2; % Bandpass-filter bandwidth

zIP_ratio_sim = zeros(N_trial, zIP_params.window_size + 1);