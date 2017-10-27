%% Depth estimaiton on synthetic dataset
%  author: Ying Jia
%  contact: jyustc93@mail.ustc.edu.cn
%  paper: ICMEW2017
%  terms of use : 
%  Any scientific work that makes use of our code should appropriately
%  mention this in the text and cite our ICMEW 2017 paper. For commercial
%  use, please contact us.
%  bibtex to cite:
%  @inproceedings{jia2017multi,
    title={Multi-occlusion handling in depth estimation of light fields},
    author={Jia, Ying and Li, Weihai},
    booktitle={Multimedia \& Expo Workshops (ICMEW), 2017 IEEE International Conference on},
    pages={13--18},
    year={2017},
    organization={IEEE}
    }
%  primary references:
%  Ting-Chun Wang, Alexei A. Efros, and Ravi Ramamoorthi.
%  Occlusion-aware depth estimation using light-field cameras. 
%  In Proceedings of International Conference on Computer Vision (ICCV), 2015.
clc
clear
cd('F:/ICME2017/Code/required/mex');
mex ./occlusionap1_mex.cpp
cd('..');
cd('..');
%% read Wanner's data 
hinfo_data = hdf5info('input/wanner_dataset/Mona.h5');
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
%  hinfo_data = hdf5info('input/wang_dataset/bedroom.h5');
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
