addpath(genpath([pwd,'\zIP']))
addpath(genpath([pwd,'\model']))
addpath(genpath([pwd,'\data_processing']))
addpath(genpath([pwd,'\Santos2016Data']))
warning('off','MATLAB:table:ModifiedAndSavedVarnames');

% Set figure properties
set(0, 'DefaultLineLineWidth', 2);
set(groot,'defaultAxesFontSize',16);
set(0,'defaultfigurecolor',[1 1 1]); % white figure background