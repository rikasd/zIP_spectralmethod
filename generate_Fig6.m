% Generate subplots in Fig. 6 of paper
%
% Rika Sugimoto-Dimitrova (rikasd@mit.edu)
% 2024-02-16

run set_default_sim_params

%% Run simulations and compute zIP for...

%% different beta
fh = figure; hold on;
colors = get(gca,'colororder');
colors_sim = ...
    [0.3 0.6 0.8;
     0.9 0.5 0.4;
     0.9 0.7 0.3];
colors_analytic = ...
    [0 0.3 0.55;
    0.64 0.08 0.18;
    0.90 0.50 0.20];
input_struct.simMotorNoiseRatio = 0.9;
alpha = 1e6;
betaVec = [0.1,0.3,2];
for iBeta = 1:length(betaVec)
    controller_params.R = alpha*diag([betaVec(iBeta),1/betaVec(iBeta)]);
    input_struct.controller_params = controller_params;    
    clear simulate_nonlinDIP;
    for iTrial = 1:N_trial
        [output_struct] = simulate_nonlinDIP(input_struct);
        [f_zIP_sim, zIP_ratio] = getZIPfromData(output_struct,zIP_params);
        zIP_ratio_sim(iTrial,:) = zIP_ratio;
    end
    figure(fh.Number);
    p_beta(iBeta) = plot(f_zIP_sim(5:end),mean(zIP_ratio_sim(:,5:end)),...
        'Linewidth',2,'Color',colors_sim(iBeta,:));
    
    input_struct.f = f_zIP_sim;
    [f_zIP_model,zIP_ratio_model] = predictZIPfromModel(input_struct);
    p_analytic(iBeta) = plot(f_zIP_model(5:end),zIP_ratio_model(5:end),'--',...
        'Linewidth',2,'Color',colors_analytic(iBeta,:));

end
figure(fh.Number);
yline(1,'k--');
legend(p_beta,{'0.1','0.3','2'})
ylim([0,2.5])
xlim([0,8])
xlabel('Frequency (Hz)');
ylabel('z_{IP}/z_{CoM}');
title('Effect of \beta on z_{IP}');

%% different motor noise ratios
fh = figure; hold on;
alpha = 1e6;
beta = 0.3;
sigmaVec = [0.5,0.9,2];
for iSigma = 1:length(sigmaVec)
    controller_params.R = alpha*diag([beta,1/beta]);
    input_struct.controller_params = controller_params;
    input_struct.motorNoiseRatio = sigmaVec(iSigma);
    clear simulate_nonlinDIP
    for iTrial = 1:N_trial
        [output_struct] = simulate_nonlinDIP(input_struct);
        [f_zIP_sim, zIP_ratio] = getZIPfromData(output_struct,zIP_params);
        zIP_ratio_sim(iTrial,:) = zIP_ratio;
    end
    figure(fh.Number);hold on;
    p_sigma(iSigma) = plot(f_zIP_sim(5:end),mean(zIP_ratio_sim(:,5:end)),...
        'Linewidth',2,'Color',colors_sim(iSigma,:));
    
    input_struct.f = f_zIP_sim;
    [f_zIP_model,zIP_ratio_model] = predictZIPfromModel(input_struct);
    p_analytic(iSigma) = plot(f_zIP_model(5:end),zIP_ratio_model(5:end),'--',...
        'Linewidth',2,'Color',colors_analytic(iSigma,:));

end
figure(fh.Number);
yline(1,'k--');
legend(p_sigma,{'0.5','0.9','2'})
ylim([0,2.5])
xlim([0,8])
xlabel('Frequency (Hz)');
ylabel('z_{IP}/z_{CoM}');
title('Effect of \sigma_r on z_{IP}');