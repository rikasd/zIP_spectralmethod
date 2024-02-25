% Generate plot in Fig. 4 of paper
%
% Rika Sugimoto-Dimitrova (rikasd@mit.edu)
% 2024-02-16

CENTI2METER = 0.01;
HEIGHT2COMZ = 0.56;

FREQSAMP_FORCEPLATE_HZ = 100; % sampling frequency for force plate data

PLANES = {'sgt','frt'};

%% Load data info
datainfo = readtable('Santos2016Data\BDSinfo.xlsx');

%% Select data file index
idatafile = 88; % Subject 8 - young healthy adult (see datainfo for details)

%% Load data
datatable = ...
    readtable(['.\Santos2016Data\' datainfo.Trial{idatafile} '.txt']);

%% Select parameters
iP = 1; % index of PLANES (select plane of interest from PLANES above)

bool_bpf = true;
bool_cpsd = false;
bool_cpsd2bpf = true;

plane = PLANES{iP};

%% Get foot-force data

dt = 1/FREQSAMP_FORCEPLATE_HZ; % sampling time step

time_s = table2array(datatable(:,'Time_s_'));
% measured forces are forces exerted on the force plates by the
% subject, and is in the opposite direction to the desired 
% ground-on-foot interaction forces, gff_N
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
com_m = ones(size(cop_m))...
    .*[0;totalHeight_m*HEIGHT2COMZ];

data_forceplate.COM = com_m; 
data_forceplate.COP = cop_m;
data_forceplate.FootForce = gff_N;
data_forceplate.sampFreq_Hz = FREQSAMP_FORCEPLATE_HZ;

%% Compute zIP

ip_params.method = 'cpsd2bpf';
ip_params.window_size = length(data_forceplate.COP);
[f_cpsd2bpf_Hz, zIP_ratio_cpsd2bpf] = getZIPfromData(data_forceplate,ip_params);

ip_params.method = 'bpf';
[f_bpf_Hz, zIP_ratio_bpf] = getZIPfromData(data_forceplate,ip_params);

%% Plot zIP

figure; hold on;
plot(f_bpf_Hz,zIP_ratio_bpf,'-','MarkerSize',15);
plot(f_cpsd2bpf_Hz,zIP_ratio_cpsd2bpf,'--','MarkerSize',15);
xlabel('Frequency (Hz)');
ylabel('z_{IP}/z_{CoM}');
xlim([0,8]);
ylim([0,2.5]);
yline(1,'k--');
legend('BPF method','Spectral method (frequency averaged)');
box off