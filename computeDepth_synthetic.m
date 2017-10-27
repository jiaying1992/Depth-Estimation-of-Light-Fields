%% Function of computing depth on synthetic dataset
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
function depth_output = computeDepth_synthetic(data)
%% INTERNAL PARAMETERS

% LF sizes                        --------------
UV_diameter         = size(data, 4)                                       ;
UV_radius           = floor(UV_diameter/2)                                ;
UV_size             = UV_diameter^2                                       ;

% Shearing                        --------------
depth_resolution        = 101                                             ;
alpha_min               = -1                                              ;
alpha_max               = 1                                               ;

% Regularize                      --------------
WS_PENALTY_W1           = 0.6                                             ;
WS_PENALTY_W2           = 0.2                                             ;
lambda_flat             = 2                                               ;
lambda_smooth           = 2                                               ;
ROBUSTIFY_SMOOTHNESS    = 1                                               ;
gradient_thres          = 1.0                                             ;
SOFTEN_EPSILON          = 1.0                                             ;
CONVERGE_FRACTION       = 1                                               ;
dilate_amount           = 4                                               ; 
addpath('required')                                                       ;
addpath('required/mex')                                                   ;
addpath('required/utilities')                                             ;
addpath('required/utilities/regularize')                                  ;
addpath('regularize')                                                     ;
addpath('regularize/GCMex')                                               ;
                                                     

x_size      = size(data, 2);                                  % spatial image height
y_size      = size(data, 1);                                  % spatial image width
LF_x_size   = (x_size) * UV_diameter;                         % total image height
LF_y_size   = (y_size) * UV_diameter;                         % total image width
LF_Remap    = reshape(permute(data, ...
               [4 1 5 2 3]), [LF_y_size LF_x_size 3]);        % the remap image
IM_Pinhole  = data(:,:,:,UV_radius+1,UV_radius+1);            % the pinhole image
LF_parameters  = struct('LF_x_size',LF_x_size,...
                             'LF_y_size',LF_y_size,...
                             'x_size',x_size,...
                             'y_size',y_size,...
                             'UV_radius',UV_radius,...
                             'UV_diameter',UV_diameter,...
                             'UV_size',UV_size,...
                             'depth_resolution',depth_resolution,...
                             'alpha_min',alpha_min,...
                             'alpha_max',alpha_max,...
                             'WS_PENALTY_W1',WS_PENALTY_W1,...
                             'WS_PENALTY_W2',WS_PENALTY_W2,...
                             'lambda_flat',lambda_flat,...
                             'lambda_smooth',lambda_smooth,...
                             'ROBUSTIFY_SMOOTHNESS',ROBUSTIFY_SMOOTHNESS,...
                             'gradient_thres',gradient_thres,...
                             'SOFTEN_EPSILON',SOFTEN_EPSILON,...
                             'CONVERGE_FRACTION',CONVERGE_FRACTION...
                             );
%% Edge detection
tic; disp('1、Edge detection')
im_edge = edge(IM_Pinhole(:,:,1), 'canny') | edge(IM_Pinhole(:,:,2), 'canny') | edge(IM_Pinhole(:,:,3), 'canny');
edge1=imdilate(im_edge,strel('disk',dilate_amount));
orient = skeletonOrientation(im_edge, [5 5]);
orient(~im_edge) = -100; 
dir = imdilate(orient, strel('disk', dilate_amount));
dir=dir*pi/180;
fprintf('Edge detecion completed in %.2f sec\n', toc);
%% Mask construction
tic; disp('2、Mask consturction')
[mask1,mask2,mask3]=maskap(IM_Pinhole,LF_parameters,edge1);
fprintf('Mask construction completed in %.2f sec\n', toc);
%% Compute Defocus and Correspondence Responses                 
fprintf('3、Computing defocus and correspondence responses\n');
tic;
[depth,depth_cost] = ...
    occlusionap1_mex(x_size, y_size, UV_diameter, LF_Remap,IM_Pinhole, alpha_min, alpha_max, depth_resolution, mask1,mask2,mask3,dir);
fprintf('Computing defocus and correspondence responses completed in %.3f seconds\n',toc);
depth=double(depth);
%% Final depth estimation
tic; disp('4. Final depth estimation')
confidence   = confCompute(depth_cost, y_size, x_size);
depth_output = DEPTH_MRF2_s(depth, confidence, im_edge, IM_Pinhole, y_size, x_size, depth_resolution) / depth_resolution;
depth_output = medfilt2(depth_output, [3 3]);
depth_output_min=min(min(depth_output));
depth_output_max=max(max(depth_output));
depth_output=(depth_output-depth_output_min)/(depth_output_max-depth_output_min);
fprintf('Final depth estimation completed in %.2f sec\n', toc);
end
