%% Run all scripts saving figures

clear
close all


%% Make sure paper parent directory is in Matlab's path

%
addpath(genpath(paper_directory()))


%% Fig. 1 - large-scale map

run(fullfile(paper_directory(), 'figures', 'makefig_01.m'))


%% Fig. 2 - small-scale map around the bump

run(fullfile(paper_directory(), 'figures', 'makefig_02.m'))


%% Fig. 3 - epsilon methods

run(fullfile(paper_directory(), 'figures', 'makefig_03.m'))


%% Fig. 4 - timeseries

run(fullfile(paper_directory(), 'figures', 'makefig_04.m'))


%% Fig. 5 - model and data

run(fullfile(paper_directory(), 'figures', 'makefig_05.m'))


%% Fig. 6 - NT1 and T1 data

run(fullfile(paper_directory(), 'figures', 'makefig_06.m'))


%% Fig. 7 - Fast-CTD and T1 data

run(fullfile(paper_directory(), 'figures', 'makefig_07.m'))


%% Fig. 8 - tidal ellipses map

run(fullfile(paper_directory(), 'figures', 'makefig_08.m'))


%% Fig. 9 - knife-edge scattering setup

% % % run(fullfile(paper_directory(), 'figures', 'makefig_09.m'))


%% Fig. 10  - scaling of epsilon with velocity

run(fullfile(paper_directory(), 'figures', 'makefig_10.m'))


%%
% To crop white margins, run on the terminal
%
% pdf-crop-margins -v -p 0 figure_file.pdf
%
%
% pdfCropMargins can be installed via:
%     pip install pdfCropMargins --upgrade

