clc
clear all

% a simply CRF for binarydenoising, with gridmodel and noisy input
% path_name = './Dataset/fakeData/';
% path_name = './Dataset/europaData/entire_log_old/res50/';
% path_name = './Dataset/europaData/gridmaps/old_features/';
path_name = './Dataset/europaData/entire_log_new/res50/';
% MCRF_binarydenoising(path_name); % when using binary segmentation with fake image(for training, testing on the very same one, just to plot the results at the end)
% MCRF_binaryEuropa_oneentirelog(path_name); % when using binary segmentation with just one Europa image (for training, testing on the very same one, just to plot the results at the end)
% MCRF_binaryEuropa(path_name);  % when using binary segmentation with two or more training examples (entire log or gridmaps)
MCRF_multiEuropa(path_name); % when using multi segmentation with two or more training examples(entire log new)

%% 
% clc
% clear all

% training CRFs on the Stanford Backgrounds Dataset
% MCRF_backgrounds;

%% 
% test for HOG feature creation
%  I=double((rgb2gray(imread('./Datasets/iccv09Data/images/0000176.jpg')))); figure(1); im(I)
%  tic, H=hog(I,8,9); toc, V=hogDraw(H,25); figure(2); im(V)