clc
clear
file_path     = 'F:/ICME2017/Code/input/tao_dataset/IMG_2.jpg';
data          = (imread(file_path)) ;    
depth_output  = computeDepth_lytro(data)  ;
figure; 
imshow(depth_output);