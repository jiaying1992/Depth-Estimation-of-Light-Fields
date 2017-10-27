%% Depth estimation on a real world image taken by Lytro camera
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
file_path     = 'F:/ICME2017/Code/input/tao_dataset/IMG_2.jpg';
data          = (imread(file_path)) ;    
depth_output  = computeDepth_lytro(data)  ;
figure; 
imshow(depth_output);
