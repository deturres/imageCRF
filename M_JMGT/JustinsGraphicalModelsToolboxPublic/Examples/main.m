clc
clear all

% a simply CRF for binarydenoising, with gridmodel and noisy input
% path_name = './Dataset/fakeData/';
% path_name = './Dataset/europaData/entire_log/res50/'
% MCRF_binarydenoising(path_name);
path_name = './Dataset/europaData/gridmaps/';
MCRF_binaryEuropa(path_name);

%% 
% clc
% clear all

% training CRFs on the Stanford Backgrounds Dataset

% MCRF_backgrounds;

%% 
% test for HOG feature creation
%  I=double((rgb2gray(imread('./Datasets/iccv09Data/images/0000176.jpg')))); figure(1); im(I)
%  tic, H=hog(I,8,9); toc, V=hogDraw(H,25); figure(2); im(V)