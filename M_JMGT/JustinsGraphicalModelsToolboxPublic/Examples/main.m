clc
clear all

% a simply CRF for binarydenoising, with gridmodel and noisy input

% OLD FEATURES
% path_name = './Dataset/fakeData/';
% path_name = './Dataset/europaData/entire_log_old/res50/';
% path_name = './Dataset/europaData/gridmaps/old_features/';

% NEW FEATURES
% path_name = './Dataset/europaData/entire_log_new/res10/';
path_name = './Dataset/europaData/entire_log_new/res10_gridmap/';

% when using BINARY segmentation with fake images(for training, testing on the very same one, just to plot the results at the end)
% MCRF_binarydenoising(path_name);

% when using BINARY segmentation with just one Europa images(for training, testing on the very same one, just to plot the results at the end)
% MCRF_binaryEuropa_entirelog_old(path_name); 

% when using BINARY segmentation with two or more training examples and different testing dataset(entire log old or gridmaps)
% [tried with new feats, really bad!)
% MCRF_binaryEuropa(path_name);  

% when using MULTI segmentation with two or more training examples and different testing dataset(entire log new[entire image or small portion/just curb&co or area(4CLASSES)])
% MCRF_multiEuropa(path_name); 

% when using MULTI segmentation with MULTI FEATURES FROM SIDEWALKDETECTOR with two or more training examples and different testing dataset(entire log new[entire image or small portion/just curb&co or area(4CLASSES)])
% MCRF_multiEuropa_multifeatures(path_name);% when using MULTI segmentation with MULTI FEATURES FROM SIDEWALKDETECTOR with two or more training examples and different testing dataset(entire log new[entire image or small portion/just curb&co or area(4CLASSES)])

% when using MULTI segmentation with MULTI FEATURES FROM SIDEWALKDETECTOR with two or more training examples and different testing dataset(entire log new[gridmap from code, import with .dat files, area(4CLASSES)])
MCRF_multiEuropa_multifeatures_gridmaps(path_name); 
%% 
% clc
% clear all

% training CRFs on the Stanford Backgrounds Dataset
% MCRF_backgrounds;

%% test for HOG feature creation
%  I=double((rgb2gray(imread('./Datasets/iccv09Data/images/0000176.jpg')))); figure(1); im(I)
%  tic, H=hog(I,8,9); toc, V=hogDraw(H,25); figure(2); im(V)

