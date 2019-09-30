%example code
clear all;

%add pyramid toolbox to path
addpath('E:\Pavan\pano_test\matlabPyrTools-master');

patch_size = 100;

%load constituent images
for k = 1:5
    im_constituent{k} = imread(['images\' num2str(k) '.jpg']);
end

%load stitched image
stitched_im = imread('images\stitched.jpg');

%% remove overlapping regions from one of the images
% This part of the code is optional. It is used to extract regions such that
% overlapping regions are used only once for processing. This can be done 
% manually or using the code below. The below code uses SIFT to remove out 
% overlapping regions. If the images have large resolution this can be 
% quite time consuming

%add vl_feat
% vl_feat can be downloaded from http://www.vlfeat.org/index.html

ip = pwd;
cd E:\Pavan\pano_test\vlfeat-0.9.19\toolbox;
feval('vl_setup');
cd(ip);
addpath('modelspecific');
addpath('multigs');

im_constituent = extract_novp(im_constituent);

%% Feature extraction
%extract features from constituent images
feat_c = extract_feat_constituent(im_constituent,patch_size);

%extract features from stitched image
feat_g = extract_feat_stitched(stitched_im,patch_size);

%% score prediction
feat = feat_c - feat_g;feat = feat';

%feature normalization
load('model_parameters');
load('model');
feat = (2./(high-low)).*(feat - (high+low)/2);

% add libsvm to path
% change path for linux and Mac
addpath('E:\Pavan\pano_test\matlabPyrTools-master\libsvm-3.22\windows');
quality_score = svmpredict(randi(100),feat,model,'-q')';