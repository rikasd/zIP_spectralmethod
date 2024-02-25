% Generate subplots in Fig. C1 of paper appendix
%
% Rika Sugimoto-Dimitrova (rikasd@mit.edu)
% 2024-02-16

iRAW = 1;
iZEROMEAN = 2;
iDETREND = 3;

iRECTWIN = 1;
iHANNWIN = 2;

WINDOWSIZES = 2.^flip(6:11);

CENTI2METER = 0.01;
HEIGHT2COMZ = 0.56;

FREQSAMP_FORCEPLATE_HZ = 100; % sampling frequency for force plate data

PLANES = {'sgt','frt'};

ip_params.method = 'cpsd';
legend_str = {'20 s       (12)','10 s       (30)','5.1 s      (66)',...
    '2.6 s      (135)','1.3 s      (276)','0.64 s    (558)'};

%% Get simulated data

% Set simulation parameters
run set_default_sim_params

N_trials = 3;

% Run simulation
dt = 1/input_struct.sampFreq_Hz;          % sampling time step
time_s = 0:dt:input_struct.simDuration_s;

clear simulate_nonlinDIP;
for iT = 1:N_trials
    [output_struct] = simulate_nonlinDIP(input_struct);
    data_SIM{iT}.COM = output_struct.COM;
    data_SIM{iT}.COP = output_struct.COP;
    data_SIM{iT}.FootForce = output_struct.FootForce;
    data_SIM{iT}.sampFreq_Hz = input_struct.sampFreq_Hz;
end % for iT
clear simulate_nonlinDIP;
        
%% Get human data

% Load data info
datainfo = readtable('Santos2016Data\BDSinfo.xlsx');

count_trial = 0;

% Select data file indices
for idatafile = [88, 89, 90] % Subject 8: young healthy female subject
                             % eyes open, firm ground 
                             % (see datainfo for details)

    % Load data
    datatable = ...
        readtable(['.\Santos2016Data\' datainfo.Trial{idatafile} '.txt']);
    % Select parameters
    iP = 1; % index of PLANES (select plane of interest from PLANES above)
    plane = PLANES{iP};
    
    % Get foot-force data
    dt = 1/FREQSAMP_FORCEPLATE_HZ; % sampling time step
    
    time_s = table2array(datatable(:,'Time_s_'));
    measuredforces_N = table2array(datatable(:,{'Fx_N_','Fy_N_','Fz_N_'}));
    cop_temp_m = table2array(datatable(:,{'COPx_cm_','COPy_cm_'}))*CENTI2METER;
    
    gff_N = zeros(2,length(measuredforces_N)); % ground-on-foot interaction forces
    cop_m = zeros(1,length(cop_temp_m));       % CoP horizontal position
    if strcmp(plane,'sgt')
        % AP direction: (+) front, (-) back
        gff_N(1,:) = (-1)*measuredforces_N(:,1);
        cop_m(1,:) = cop_temp_m(:,1);
    elseif strcmp(plane, 'frt')
        % ML direction: (+) left, (-) right
        gff_N(1,:) = (-1)*measuredforces_N(:,2);
        cop_m(1,:) = cop_temp_m(:,2);
    end        
    gff_N(2,:) = measuredforces_N(:,3); % (+) up, (-) down
    
    totalHeight_m = datainfo.Height(idatafile)*CENTI2METER;
    com_m = ones(size(cop_m)).*[0;totalHeight_m*HEIGHT2COMZ];
    
    count_trial = count_trial + 1;
    data_Santos{count_trial}.COM = com_m; 
    data_Santos{count_trial}.COP = cop_m;
    data_Santos{count_trial}.FootForce = gff_N;
    data_Santos{count_trial}.sampFreq_Hz = FREQSAMP_FORCEPLATE_HZ;

end % for idatafile

%% Plot zIP curve
clear zIP_ratio_cpsd f_cpsd_Hz

data_All = {data_SIM, data_Santos};

for idata = 1:2

    indices_winsize = 1:length(WINDOWSIZES);
    indices_detrend = iDETREND;%[iRAW,iZEROMEAN,iDETREND];
    indices_winfilt = iHANNWIN;%[iRECTWIN,iHANNWIN];
    
    index_init = 1; % skip the first N data points
    
    zIP_ratio_cpsd = ...
        cell(length(indices_winsize),length(indices_detrend),length(indices_winfilt));
    f_cpsd_Hz = zIP_ratio_cpsd;
    
    for iWinsize = indices_winsize
    for iDetrend = indices_detrend
    for iWinfilt = indices_winfilt
    
    ip_params.window_size = WINDOWSIZES(iWinsize);
    ip_params.detrend_option = iDetrend;
    ip_params.window_taper_option = iWinfilt;
    
    [f_zIP_Hz, zIP_ratio_cpsd_temp] = getZIPfromData(data_All{idata},ip_params);
    zIP_ratio_cpsd{iWinsize,iDetrend,iWinfilt} = zIP_ratio_cpsd_temp;
    f_cpsd_Hz{iWinsize,iDetrend,iWinfilt} = f_zIP_Hz;
    
    end % iWinfilt
    end % iDetrend
    end % iWinsize
    
    fh = figure;hold on;
    colors = get(gca,'colororder');    
    lineSpecs.width = 1.2;
    lineSpecs.alpha = 0.5;
    lineSpecs.style = '-';
    icolor = 1;
    for iWinsize = indices_winsize
    for iDetrend = indices_detrend
    for iWinfilt = indices_winfilt
        lineSpecs.color = colors(icolor,:);
        plot(f_cpsd_Hz{iWinsize,iDetrend,iWinfilt}(index_init:end),...
            zIP_ratio_cpsd{iWinsize,iDetrend,iWinfilt}(index_init:end),...
            'Color',lineSpecs.color,'Linewidth',lineSpecs.width);
        icolor = icolor + 1;
    end % iWinfilt
    end % iDetrend
    end % iWinsize
    yline(1,'k--');
    xlabel('Frequency (Hz)');xlim([0,8]);
    ylabel('z_{IP}/z_{CoM}');ylim([0,2.5]);
    set(gca,'Fontsize',15);
    legend(legend_str);
    box off

end % idata