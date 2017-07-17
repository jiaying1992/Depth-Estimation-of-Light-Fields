#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <algorithm>
#include "mex.h"
using namespace std;

#define Depth        plhs[0]
#define Depth_cost   plhs[1]



#define X_size_fin         prhs[0]
#define Y_size_fin         prhs[1]
#define UV_diameter        prhs[2]
#define Im_in_remap        prhs[3]
#define Im_pinhole         prhs[4]
#define Alpha_min          prhs[5]
#define Alpha_max          prhs[6]
#define Depth_res          prhs[7]
#define Mask_ap1           prhs[8]
#define Mask_ap2           prhs[9]
#define Mask_ap3           prhs[10]
#define Orient             prhs[11]


double alpha_min;
double alpha_max;
int depth_res;
double alpha_step;
int width;
int height;
int pixelNum;
int uv_diameter;
int uv_radius;
int uv_size;
int remap_height;
int remap_width;
int remap_pixelNum;
double ratio = 1;          // weight of refocus cue to correspondence cue
int bsz = 5;               // spatial window size of bilateral filter for depth cost
double sigma_range = 0.01; // range sigma of bilateral filter
double sigma_c = 0.01;
double sigma_f = 0.01;

int checkX(int x) {
    if (x < 0) return 0;
    else if (x >= width) return width-1;
    else return x;
}

int checkY(int y) {
    if (y < 0) return 0;
    else if (y >= height) return height-1;
    else return y;
}

double getVar(double array[], int size, double* mean)
{
    *mean = 0;
    int num = 0;
    for (int i = 0; i < size; i++) {
        if (array[i] >= 0) {
            *mean += array[i];
            num++;
        }
    }
    *mean /= num;
    double var = 0;
    for (int i = 0; i < size; i++) {
        if (array[i] >= 0) {
            var += (array[i]-*mean) * (array[i]-*mean);
        }
    }
    if (num == 1) var = INT_MAX;
    else var /= (num-1);
    return var;
}
double getVar1(double array[],int size,double im)
{
      int num=0;
      double var1=0;
      for(int i=0;i<size;i++){
         if(array[i]>=0){
             var1+=(array[i]-im)*(array[i]-im);
             num++;
          }
       }
       if(num==1) var1=INT_MAX;
       else var1 /= (num-1);
       return var1;
}



void remapping(double im_in_remap[], int x, int y, int alpha_num, double remap[])
{
    double alpha = alpha_min + alpha_step * alpha_num;
    for (int i = -uv_radius; i <= uv_radius; i++) {
        for (int j = -uv_radius; j <= uv_radius; j++) {
            float y_ind   = (float)i*alpha + y;
            float x_ind   = (float)j*alpha + x;
            
            int x_floor = floor(x_ind);
            int y_floor = floor(y_ind);
            
            int x_1     = checkX(x_floor);
            int y_1     = checkY(y_floor);
            int x_2     = checkX(x_floor+1);
            int y_2     = checkY(y_floor+1);
            
            float x_1_w   = 1-(x_ind-x_floor);
            float x_2_w   = 1-x_1_w;
            float y_1_w   = 1-(y_ind-y_floor);
            float y_2_w   = 1-y_1_w;
            
            int x_1_index = j+uv_radius + (x_1)*uv_diameter;
            int y_1_index = i+uv_radius + (y_1)*uv_diameter;
            int x_2_index = j+uv_radius + (x_2)*uv_diameter;
            int y_2_index = i+uv_radius + (y_2)*uv_diameter;
            
            for (int c = 0; c < 3; c++) {
                remap[(i+uv_radius)*uv_diameter+(j+uv_radius)+c*uv_size] =
                y_1_w * x_1_w * im_in_remap[y_1_index+x_1_index*remap_height+c*remap_pixelNum] +
                y_2_w * x_1_w * im_in_remap[y_2_index+x_1_index*remap_height+c*remap_pixelNum] +
                y_1_w * x_2_w * im_in_remap[y_1_index+x_2_index*remap_height+c*remap_pixelNum] +
                y_2_w * x_2_w * im_in_remap[y_2_index+x_2_index*remap_height+c*remap_pixelNum];
            }
        }
    }
}

void remap_kmeans(double im_in_remap[],double im_pinhole[],double mask_ap1[],double mask_ap2[],int x,int y,float defocus_response[],float corresp_response[],double orient[],double mask_ap3[])
{
	double* remap    = new double[uv_size*3];
    double* remap1   = new double[uv_size*3];
    double* remap2   = new double[uv_size*3];
    double* remap3   = new double[uv_size*3];
    double* remap_p1 = new double[uv_size*3];
    double* remap_p2 = new double[uv_size*3];
	double mean_p1[3];
    double mean_p2[3];
    double mean_p3[3];
    double mean_p4[3];
    double mean_p5[3];
	for (int alpha_num=0;alpha_num<depth_res;alpha_num++){
		remapping(im_in_remap,x,y,alpha_num,remap);
		for(int c=0;c<3;c++){
			for(int i=0;i<uv_diameter;i++){
				for(int j=0;j<uv_diameter;j++){
					remap1[j+i*uv_diameter+c*uv_size]=(mask_ap1[x*remap_height*uv_diameter+i*remap_height+y*uv_diameter+j]==1) ? remap[j+i*uv_diameter+c*uv_size]:-1;
                    remap2[j+i*uv_diameter+c*uv_size]=(mask_ap2[x*remap_height*uv_diameter+i*remap_height+y*uv_diameter+j]==1) ? remap[j+i*uv_diameter+c*uv_size]:-1;
                    remap3[j+i*uv_diameter+c*uv_size]=(mask_ap3[x*remap_height*uv_diameter+i*remap_height+y*uv_diameter+j]==1) ? remap[j+i*uv_diameter+c*uv_size]:-1;
				}
			}
		}
        double var1 = getVar(remap1, uv_size, mean_p1)+getVar(&remap1[uv_size], uv_size, &mean_p1[1])+getVar(&remap1[uv_size*2], uv_size, &mean_p1[2]);
        double var2 = getVar1(remap1, uv_size, im_pinhole[y+x*height])+ getVar1(&remap1[uv_size], uv_size, im_pinhole[y+x*height+pixelNum])+getVar1(&remap1[uv_size*2], uv_size, im_pinhole[y+x*height+2*pixelNum]); 
        double var3 = getVar(remap2, uv_size, mean_p2)+getVar(&remap2[uv_size], uv_size, &mean_p2[1])+getVar(&remap2[uv_size*2], uv_size, &mean_p2[2]);
        double var4 = getVar1(remap2, uv_size, im_pinhole[y+x*height])+ getVar1(&remap2[uv_size], uv_size, im_pinhole[y+x*height+pixelNum])+getVar1(&remap2[uv_size*2], uv_size, im_pinhole[y+x*height+2*pixelNum]);
        double var5 = getVar(remap3, uv_size, mean_p3) + getVar(&remap3[uv_size], uv_size, &mean_p3[1]) + getVar(&remap3[uv_size*2], uv_size, &mean_p3[2]);
        double var6 = getVar1(remap3, uv_size, im_pinhole[y+x*height])+ getVar1(&remap3[uv_size], uv_size, im_pinhole[y+x*height+pixelNum])+getVar1(&remap3[uv_size*2], uv_size, im_pinhole[y+x*height+2*pixelNum]);
        double var7=min(sqrt(var1),sqrt(var3));
        double var8=min(sqrt(var2),sqrt(var4));
        double var9=min(var7,sqrt(var5));
        double var10=min(var8,sqrt(var6));
        
        if(orient[y+x*height] >=-100*M_PI/180+0.1){
            double theta = orient[y+x*height];
            for (int c = 0; c < 3; c++) {
               for (int i = 0; i < uv_diameter; i++) {
                  for (int j = 0; j < uv_diameter; j++) {
                    remap_p1[i*uv_diameter + j + c*uv_size] = (uv_radius-i-sin(theta+M_PI/2) > (tan(theta) * (j-cos(theta+M_PI/2)-uv_radius)))?
                                                                remap[i*uv_diameter + j + c*uv_size] : -1;
                    remap_p2[i*uv_diameter + j + c*uv_size] = (uv_radius-i-sin(theta+M_PI/2) < (tan(theta) * (j-cos(theta+M_PI/2)-uv_radius)))?
                                                                remap[i*uv_diameter + j + c*uv_size] : -1;
                  }
                }
            }
        // compute both the correspondence cue: patch variance, and refocus cue: patch mean
        double var11 = getVar(remap_p1, uv_size, mean_p4) + getVar(&remap_p1[uv_size], uv_size, &mean_p4[1]) + getVar(&remap_p2[uv_size*2], uv_size, &mean_p4[2]);
        double var12 = getVar1(remap_p1, uv_size, im_pinhole[y+x*height])+ getVar1(&remap_p1[uv_size], uv_size, im_pinhole[y+x*height+pixelNum])+getVar1(&remap_p1[uv_size*2], uv_size, im_pinhole[y+x*height+2*pixelNum]); 
        double var13 = getVar(remap_p2, uv_size, mean_p5) + getVar(&remap_p2[uv_size], uv_size, &mean_p5[1]) + getVar(&remap_p2[uv_size*2], uv_size, &mean_p5[2]);
        double var14 = getVar1(remap_p2, uv_size, im_pinhole[y+x*height])+ getVar1(&remap_p2[uv_size], uv_size, im_pinhole[y+x*height+pixelNum])+getVar1(&remap_p2[uv_size*2], uv_size, im_pinhole[y+x*height+2*pixelNum]);
        double var15 = min(sqrt(var11),sqrt(var13));
        double var16 = min(sqrt(var12),sqrt(var14));
        defocus_response[y + x*height + alpha_num*pixelNum]= min(var10,var16);
        corresp_response[y + x*height + alpha_num*pixelNum] = min(var9,var15);
       }
        else
        {
            defocus_response[y + x*height + alpha_num*pixelNum]=  var10;
            corresp_response[y + x*height + alpha_num*pixelNum] = var9;
        }
        
	}
	delete[] remap;
    delete[] remap1;
    delete[] remap2;
    delete[] remap3;
    delete[] remap_p1;
    delete[] remap_p2;
   
}

void costFilter( float corresp_response[],float defocus_response[], float corresp_cost[],float defocus_cost[], double im_pinhole[])
{
    for (int i = 0; i < depth_res; i++) {
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                double sumVal_c = 0, sumVal_f = 0;
                double weight_sum_c = 0, weight_sum_f = 0;
                for (int u = -bsz/2; u <= bsz/2; u++) {
                    for (int v = -bsz/2; v <= bsz/2; v++) {
                        int xu = checkX(x+u);
                        int yv = checkY(y+v);
                        
                        // compute bilateral weight
                        double color_dif = 0;
                        for (int c = 0; c < 3; c++)
                            color_dif = pow(im_pinhole[y+x*height+c*pixelNum]-im_pinhole[yv+xu*height+c*pixelNum], 2);
                        double weight = exp(-color_dif / (2*pow(sigma_range, 2)));
                        
                        // correspondence cost
                        sumVal_c += weight * corresp_response[xu*height+yv+i*pixelNum];
                        weight_sum_c += weight;
                        
                        // refocus cost
                        sumVal_f += weight * defocus_response[xu*height+yv+i*pixelNum];
                        weight_sum_f += weight;
                    }
                }
                sumVal_c = weight_sum_c != 0? sumVal_c / weight_sum_c : 1;
                sumVal_f = weight_sum_f != 0? sumVal_f / weight_sum_f : 1;
                corresp_cost[y+x*height+i*pixelNum]=sumVal_c;
                defocus_cost[y+x*height+i*pixelNum]=sumVal_f;
                
            }
        }
    }
}

void computequan( float depth_cost[],float corresp_cost[],float defocus_cost[])
{
    for(int y=0;y<height;y++){
        for(int x=0;x<width;x++){
            double minval_c=corresp_cost[y+x*height];
            double minval_f=defocus_cost[y+x*height];
            for(int i=1;i<depth_res;i++){
                if(corresp_cost[y+x*height+i*pixelNum]<minval_c){
                    minval_c=corresp_cost[y+x*height+i*pixelNum];
                }
                if(defocus_cost[y+x*height+i*pixelNum]<minval_f){
                    minval_f=defocus_cost[y+x*height+i*pixelNum];
                }
            }
            double quan_c=0;
            double quan_f=0;
            for(int i=0;i<depth_res;i++){
                quan_c+=exp(-pow((corresp_cost[y+x*height+i*pixelNum]-minval_c),2) / (2*pow(sigma_c, 2)));
                quan_f+=exp(-pow((defocus_cost[y+x*height+i*pixelNum]-minval_f),2) / (2*pow(sigma_f, 2)));
            }
            quan_c=1/quan_c;
            quan_f=1/quan_f;
            for(int i=0;i<depth_res;i++){
                depth_cost[y+x*height+i*pixelNum] = (quan_c/(quan_c+quan_f))*corresp_cost[y+x*height+i*pixelNum] + (quan_f/(quan_c+quan_f))*defocus_cost[y+x*height+i*pixelNum];
            }
            
        }
    }
}


        
int getMinIdx(int depth[], float depth_cost[],  int x, int y)
{
    int minIdx;
    double minVal = INT_MAX;
    for (int i = 0; i < depth_res; i++) {
        double sumVal = depth_cost[y+x*height+i*pixelNum];
        if (sumVal < minVal) {
            minVal = sumVal;
            minIdx = i;
        }
    }
    return minIdx;
}
void computeResponse(double* im_in_remap,double* mask_ap1,double* mask_ap2,double* im_pinhole,int* depth,float* depth_cost,double* orient,double* mask_ap3)
{
    float* defocus_response= new float [pixelNum*depth_res];
    float* corresp_response= new float [pixelNum*depth_res];
    float* defocus_cost=new float[pixelNum*depth_res];
    float* corresp_cost=new float[pixelNum*depth_res];
	for(int y=0;y<height;y++){
		for(int x=0;x<width;x++){
			remap_kmeans(im_in_remap,im_pinhole,mask_ap1,mask_ap2,x,y,defocus_response,corresp_response,orient,mask_ap3);
		}
	}
    // apply a bilateral filter on the depth cost
    costFilter( corresp_response,defocus_response, corresp_cost,defocus_cost, im_pinhole);
     computequan(depth_cost,corresp_cost,defocus_cost);
    // get the depth with the lowest depth cost
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            depth[x*height+y] = getMinIdx(depth, depth_cost, x, y);
        }
    }
    delete[] defocus_response;
    delete[] corresp_response;
    delete[] defocus_cost;
    delete[] corresp_cost;
}

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    double* x_size_fin_pt   = mxGetPr(X_size_fin);
	double* y_size_fin_pt   = mxGetPr(Y_size_fin);
	double* uv_diameter_pt  = mxGetPr(UV_diameter);
    double* im_in_remap_pt  = mxGetPr(Im_in_remap);
    double* im_pinhole_pt   = mxGetPr(Im_pinhole);
    double* alpha_min_pt    = mxGetPr(Alpha_min);
    double* alpha_max_pt    = mxGetPr(Alpha_max);
    double* depth_res_pt    = mxGetPr(Depth_res);
	double* mask_ap1_pt     = mxGetPr(Mask_ap1);
    double* mask_ap2_pt     = mxGetPr(Mask_ap2);
    double* mask_ap3_pt     = mxGetPr(Mask_ap3);
    double* orient_pt       = mxGetPr(Orient);
   
    int* depth_pt;
    float* depth_cost_pt;
	if(nrhs!=12)
		mexErrMsgTxt("Twelve input arguments required.");
	else if(nlhs>2)
		mexErrMsgTxt("Too many output arguments.");
	
	alpha_min = *alpha_min_pt;
	alpha_max = *alpha_max_pt;
	depth_res = *depth_res_pt;
	alpha_step = (alpha_max-alpha_min) / (depth_res-1);
	
	uv_diameter = *uv_diameter_pt;
	uv_radius = uv_diameter/2;
	uv_size = uv_diameter * uv_diameter;
	width = *x_size_fin_pt;
	height = *y_size_fin_pt;
	pixelNum = width*height;
	remap_height = height*uv_diameter;
	remap_width = width*uv_diameter;
	remap_pixelNum = remap_height*remap_width;
    
    Depth = mxCreateNumericMatrix(height, width, mxINT32_CLASS, mxREAL);
    depth_pt = (int*)mxGetData(Depth);
	
    const mwSize dims2[] = { height, width, depth_res };
    Depth_cost= mxCreateNumericArray(3, dims2, mxSINGLE_CLASS, mxREAL);
    depth_cost_pt = (float*)mxGetPr(Depth_cost);
    
	
	
	computeResponse(im_in_remap_pt,mask_ap1_pt,mask_ap2_pt,im_pinhole_pt,depth_pt,depth_cost_pt,orient_pt,mask_ap3_pt);
	
	
	}