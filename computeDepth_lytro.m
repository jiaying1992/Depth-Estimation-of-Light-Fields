%% Function of computing depth from a real world image
%  author: Ying Jia
%  contact: jyustc93@mail.ustc.edu.cn
%  time: 20170718
%  paper: ICMEW2017
%  TERMS OF USE : 
%  Any scientific work that makes use of our code should appropriately
%  mention this in the text and cite our ICMEW 2017 paper. For commercial
%  use, please contact us.
function depth_output = computeDepth_lytro(data)
addpath('required')                                                       ;
addpath('required/mex')                                                   ;
addpath('required/utilities')                                             ;
addpath('required/utilities/regularize')                                  ;
addpath('regularize')                                                     ;
addpath('regularize/GCMex')                                               ;
load('required/camera_data/image_cords')                                  ;
%% Parameters
alpha_max = 1;          
alpha_min=-1;
depth_resolution = 101;        % depth resolution
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

% parameters from input
UV_diameter = 7;                           % angular resolution
UV_radius   = 3;                    % half angular resolution
UV_size=49;
y_size           = 362;                           % spatial image height
x_size           = 311;                           % spatial image width
LF_y_size   = y_size * UV_diameter;                         % total image height
LF_x_size   = x_size * UV_diameter;                         % total image width
LF_parameters       = struct('LF_x_size',LF_x_size,...
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
                             )                                            ;

Lytro_RAW_Demosaic  = im2double(demosaic(data,'bggr'))  ;   

%% JPEG -> Our LF_Remap Standard and Pinhole Image
fprintf('1. Remapping LF JPEG to our standard                  *******\n');
tic                                                                       ;
% RAW to Remap
LF_Remap            = RAW2REMAP...
                        (Lytro_RAW_Demosaic,image_cords,LF_parameters)    ;
% Remape to Pinhole (Center View)
IM_Pinhole          = REMAP2PINHOLE...
                        (LF_Remap,LF_parameters)                          ;
fprintf('Remapping LF JPEG to our standard completed in %.3f seconds\n',toc);
%% Edge detection and orientation computation
tic; disp('2. Edge detection')
im_edge = edge(IM_Pinhole(:,:,1), 'canny') | edge(IM_Pinhole(:,:,2), 'canny') | edge(IM_Pinhole(:,:,3), 'canny');
edge1=imdilate(im_edge,strel('disk',dilate_amount));
orient = skeletonOrientation(im_edge, [5 5]);
orient(~im_edge) = -100; 
dir = imdilate(orient, strel('disk', dilate_amount));
dir=dir*pi/180;
fprintf('Edge detecion completed in %.2f sec\n', toc);
%% Mask construction
 tic; disp('3. Mask construction')
[mask1,mask2,mask3]=maskap(IM_Pinhole,LF_parameters,edge1); 
fprintf('Mask construction completed in %.2f sec\n', toc);
%% II. Compute Defocus and Correspondence Responses                 
fprintf('4. Computing Defocus and Correspondence Responses    *******\n');
tic                                                                    
[depth,depth_cost] = ...
    occlusionap1_mex(x_size, y_size, UV_diameter, LF_Remap,IM_Pinhole, alpha_min, alpha_max, depth_resolution, mask1,mask2,mask3,dir);
fprintf('Computing Defocus and Correspondence Responses completed in %.3f seconds\n',toc);
depth=double(depth);
%% Final depth estimation
tic; disp('5. Final depth estimation')
confidence   = confCompute(depth_cost, y_size, x_size);
depth_output = DEPTH_MRF2_s(depth, confidence, im_edge, IM_Pinhole, y_size, x_size, depth_resolution) / depth_resolution;
depth_output = medfilt2(depth_output, [3 3]);
depth_output_min=min(min(depth_output));
depth_output_max=max(max(depth_output));
depth_output=(depth_output-depth_output_min)/(depth_output_max-depth_output_min);
fprintf('Final depth estimation completed in %.2f sec\n', toc);
end
