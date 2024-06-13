% Main script
%
% INSTRUCTIONS:
% 1. Download the public data set from Santos and Duarte (2016)
%    and add it to the "Santos2016Data" folder.
% 2. Run this script.

clear; clc; close all;

run setup

%% Run simulations and compute zIP to generate figures from paper:

run generate_Fig5
run generate_Fig6  %<< May take some time running simulations
run generate_FigC1