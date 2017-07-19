%% Depth estimation on a real world image taken by Lytro camera
%  author: Ying Jia
%  contact: jyustc93@mail.ustc.edu.cn
%  time: 20170718
%  paper: ICMEW2017
%  TERMS OF USE : 
%  Any scientific work that makes use of our code should appropriately
%  mention this in the text and cite our ICMEW 2017 paper. For commercial
%  use, please contact us.
clc
clear
file_path     = 'F:/ICME2017/Code/input/tao_dataset/IMG_2.jpg';
data          = (imread(file_path)) ;    
depth_output  = computeDepth_lytro(data)  ;
figure; 
imshow(depth_output);
