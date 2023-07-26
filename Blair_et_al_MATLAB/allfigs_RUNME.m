% step1 - read the README.md !
clear;
% actual top directory where _DATA and _MATLAB are found
repo_dir = 'D:\Sample Data\Linear Maze Data\tadblair-main\';

addpath(repo_dir)
% folder where data files are found
topdir = [repo_dir 'Blair_et_al_DATA\'];

% create figure 1
Figure1_behavior;

% create figure 2 data
Figure2_analysis;
% create figure 2
Figure2_plotting;

% create figure 3 data
Figure3_analysis;
% create figure 4
Figure3_plotting;

% prepare figure 4 data
Figure4A_analysis;
% create figure 4
Figure4A_plotting;

% if problems or issues, check the README.md first, then contact tadblair@ucla.edu or gblair@ucla.edu