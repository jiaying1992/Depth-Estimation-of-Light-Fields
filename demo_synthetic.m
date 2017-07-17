clc
clear

cd('F:/ICME2017/Code/required/mex');
mex ./occlusionap1_mex.cpp
cd('..');
cd('..');
%% read Wanner's data 
hinfo_data = hdf5info('input/wanner_dataset/Buddha.h5');
data = hdf5read(hinfo_data.GroupHierarchy.Datasets(2));
data = permute(data, [3 2 1 5 4] );   
data = im2double(data(:, :, :, :, end:-1:1));
gt = hdf5read(hinfo_data.GroupHierarchy.Datasets(1));
gtdepth = gt(:,:,5,5);
gtdepth=im2double(gtdepth');
gtmin=min(min(gtdepth));
gtmax=max(max(gtdepth));
gtdepth=(gtdepth-gtmin)/(gtmax-gtmin);
%% read Wang's data
%  hinfo_data = hdf5info('input/wang_dataset/sculptures.h5');
%  data = double(hdf5read(hinfo_data.GroupHierarchy.Datasets(3)));
%  gt = hdf5read(hinfo_data.GroupHierarchy.Datasets(1));
%  gtdepth = gt(:,:,5,5);
%  gtdepth=im2double(gtdepth);
%  gtmin=min(min(gtdepth));  
%  gtmax=max(max(gtdepth));
%  gtdepth=(gtdepth-gtmin)/(gtmax-gtmin);
%  te=ones(600,800);
%  gtdepth=te-gtdepth;
depth_output = computeDepth_synthetic(data);
%% compute MSE
a=(depth_output-gtdepth).*(depth_output-gtdepth);
b=sum(sum(a));
