/*
    CS 89/189 Computational Aspects of Digital Photography C++ basecode.

    Adapted from MIT's 6.815/6.865 basecode, written and designed by:
        Frédo Durand
        Katherine L. Bouman
        Gaurav Chaurasia
        Adrian Vasile Dalca
        Neal Wadhwa

    With additions & modifications by
        Wojciech Jarosz
    for Dartmouth's CS 89/189.
*/

// filtering.cpp
// Assignment 4

#include "filtering.h"
#include "basicImageManipulation.h"
#include <math.h>
#include <algorithm>

using namespace std;

/**************************************************************
 //            IMAGE CONVOLUTION AND FILTERING               //
 *************************************************************/


// convolve an image with a box filter of size k by k
FloatImage boxBlur(const FloatImage &im, const int &k, bool clamp) {
    
    // create a new empty image
    if ( float(k / 2) == 0 )
        throw Kinvalid();

    FloatImage imFilter(im.width(), im.height(), im.channels());
    float itvl = 0; 
    // convolve the image with the box filter
    // cout<<imFilter.width()<<"  "<< imFilter.height()<<endl;
    for (int i = 0; i < imFilter.width(); ++i)
        for (int j = 0; j < imFilter.height(); ++j)
            for (int p = 0; p < imFilter.channels(); ++p)
            {
                itvl = 0;
                for (int f = -k/2; f <= k/2; f++)
                    for (int m = -k/2; m <= k/2; m++){
                        itvl +=  im.smartAccessor(i + f,j + m,p,clamp);      
                    }
                imFilter(i,j,p) = float( itvl / (k *k));  
            }

    return imFilter;
}


// reimeplement the box filter using the filter class.
// check that your results math those in the previous function "boxBlur"
FloatImage boxBlur_filterClass(const FloatImage &im, const int &k, bool clamp) {
// use Convolve() to apply convolution
    FloatImage imgbxft(im);
    Filter box(k,k);
    for (int i=0; i < box.width * box.height; ++i){
        box.kernel[i] = float(1.0 / (k * k));
    }
    
    imgbxft = box.Convolve(imgbxft,clamp);
    return imgbxft;// CHANGEME
}


// uses a Sobel kernel to compute the horizontal and vertical
// components of the gradient of an image and returns the gradient magnitude.
FloatImage gradientMagnitude(const FloatImage &im, bool clamp){
    
    FloatImage grd(im);
    FloatImage im_h(im.width(), im.height(), im.channels());
    FloatImage im_v(im.width(), im.height(), im.channels());
    static const int arr1[] = {-1,0,1,-2,0,2,-1,0,1}; 
    static const int arr2[] = {-1,-2,-1,0,0,0,1,2,1};
    vector<float> vect_h (arr1, arr1 + sizeof(arr1) / sizeof(arr1[0]));
    vector<float> vect_v (arr2, arr2 + sizeof(arr2) / sizeof(arr2[0]));

    Filter horiz(vect_h,3,3);
    Filter vertic(vect_v,3,3); 
    // sobel filtering in x direction
    im_h = horiz.Convolve(grd,true);
    // sobel filtering in y direction
    im_v = vertic.Convolve(grd,true);
    // compute squared magnitude
    for (int i = 0; i < im.width(); i++)
        for (int j = 0; j < im.height(); j++)
            for (int k = 0; k < im.channels(); k++)
                grd(i,j,k) = sqrt(pow(im_h(i,j,k),2) + pow(im_v(i,j,k),2));
    return grd;
}

// create a vector containing the normalized values in a 1D Gaussian filter
vector<float> gauss1DFilterValues(float sigma, float truncate){
    
    // calculate the size of the filter
    // compute the un-normalized value of the gaussian
    // normalize
    int size = 1+2 * ceil(sigma * truncate);
    std::vector<float> g_kernal(size,0);
    // cout<<"fine = "<<g_kernal.size()<<"sz = " <<size<<endl;
    float rslt = 0;
    float sum = 0; 
    for (int i = - ceil(sigma * truncate); i <= ceil(sigma * truncate); i++){  //TODO: IS THIS SAMPLE METHOD RIGHT ? 
        rslt = exp(- (i * i) / (2 * sigma * sigma));
        g_kernal[i+ceil(sigma * truncate)] = rslt;
        sum += rslt;
        // cout << "rslt = " <<rslt<<endl;
         // cout<<"kernal val = "<<g_kernal[i+2*ceil(sigma * truncate)]<<endl;
    }
    for (int j = 0; j < size; ++j)
        g_kernal[j] = g_kernal[j] / sum;

    return g_kernal;
}

// blur across the rows of an image
FloatImage gaussianBlur_horizontal(const FloatImage &im, float sigma, float truncate, bool clamp){
    FloatImage output(im);
    int height = 1;
    std::vector<float> gs_1d; 
    gs_1d = gauss1DFilterValues(sigma,truncate);
    Filter gs_h(gs_1d,gs_1d.size(),height);
    // filter in the x direction
    output = gs_h.Convolve(output,clamp);
    return output;// CHANGEME
}

// create a vector containing the normalized values in a 2D Gaussian filter
vector<float> gauss2DFilterValues(float sigma, float truncate){ // TODO: output wrong
    
    // compute the filter size
    int width = 1 + 2 * ceil(sigma * truncate);
    int height = 1 + 2 * ceil(sigma * truncate);
    int center = ceil(sigma * truncate);
    float dis = 0;
    float rslt;
    std::vector<float> gaussian_val(width * height,0);
    FloatImage gaus_fun(width,height,1);


    for (int i = 0; i < width; ++i)
        for (int j = 0; j < height; ++j)
        {
           dis = sqrt(pow(i - center,2) + pow(j - center,2)); 
            //dis = max(abs(i - center), abs(j-center)); 
            rslt = exp(- (dis * dis) / (2 * sigma * sigma)) / (sigma * sigma * 2 * M_PI);
            gaussian_val[i + j * width] = rslt;
            gaus_fun(i,j,0) = rslt;
        }

    // for (int i =0; i< width * height -1; i++)
    //     cout<< gaussian_val[i] <<endl;
    gaus_fun.write("./Output/gaus_2d.png");

    // compute the unnormalized value of the gaussian and put it in a row-major vector
    // normalize
    return gaussian_val;// CHANGEME
}


// Blur an image with a full  full 2D rotationally symmetric Gaussian kernel
FloatImage gaussianBlur_2D(const FloatImage &im, float sigma, float truncate, bool clamp){
    
    // blur using a 2D gaussian filter (use gauss2DFilterValues())
    FloatImage output(im);
    std::vector<float> g_2d = gauss2DFilterValues(sigma,truncate);
    Filter gabl_2d(g_2d,1 + 2 * ceil(sigma * truncate),1 + 2 * ceil(sigma * truncate));
    output = gabl_2d.Convolve(output,clamp);
    return output;//  CHANGEME
}

// Use principles of seperabiltity to blur an image using 2 1D Gaussian Filters
FloatImage gaussianBlur_seperable(const FloatImage &im, float sigma, float truncate, bool clamp){
    
    // blur using 2, 1D filters in the x and y directions
    FloatImage output(im);
    std::vector<float> v_val;
    
    v_val = gauss1DFilterValues(sigma,truncate);
    int width = v_val.size();
    int height = width;
    // cout<< "vector size = "<<width<<endl;
    Filter gaus_h(v_val,width,1);
    Filter gaus_v(v_val,1,height);

    output = gaus_h.Convolve(output,clamp);
    output = gaus_v.Convolve(output,clamp);


    return output;//  CHANGEME
}


// sharpen an image
FloatImage unsharpMask(const FloatImage &im, float sigma, float truncate, float strength, bool clamp){
    
    // get the low pass image and subtract it from the original image to get the high pass image
    FloatImage output(im);
    FloatImage lowpass(im.width(), im.height(), im.channels());
    FloatImage highpass(im.width(), im.height(), im.channels());

    lowpass = gaussianBlur_seperable(output,sigma,truncate,clamp);
    highpass = output - lowpass;

    output = output + strength * highpass;
    return output;
    
}


// Denoise an image using bilateral filtering
FloatImage bilateral(const FloatImage &im, float sigmaRange, float sigmaDomain, float truncateDomain, bool clamp){
    
    // calculate the filter size
    // for every pixel (x,y) in the image set value to weighted sum of values in the filter region
    int size = 1 + 2 * ceil(sigmaDomain * truncateDomain);
    float gaussian_d;
    float gaussian_r;
    float diff; 
    float overall_val; 
    float k_normalize;

    FloatImage output(im.width(),im.height(),im.channels());
    FloatImage output_inter(im);


    for (int i = 0; i < output.width(); ++i)
        for (int j = 0; j < output.height(); ++j)
            for (int k = 0; k < output.channels(); ++k)
            {
                gaussian_d = 0;
                gaussian_r = 0;
                overall_val = 0;
                k_normalize = 0;
                for (int f = - size / 2; f <= size /2; ++f)
                    for (int m = - size /2; m <= size /2; ++m){
                        gaussian_d = exp(- ((f * f) + (m * m)) / (2 * sigmaDomain * sigmaDomain)) / (sigmaDomain * sigmaDomain * 2 * M_PI);
                        diff = output_inter.smartAccessor(i,j,k) - output_inter.smartAccessor(i + f,j + m,k);
                        gaussian_r = exp(- (diff * diff) / (2 * sigmaRange * sigmaRange)) / (sigmaRange * sqrt(2 * M_PI));
                        overall_val += gaussian_d * gaussian_r * output_inter.smartAccessor(i + f,j + m,k);  // No need for flipping, Gaussian is symmetric.
                        k_normalize += gaussian_d * gaussian_r;
                    }
                output(i,j,k) = overall_val / k_normalize;
            }

    return output;//  CHANGEME
}


// Bilaterial Filter an image seperatly for
// the Y and UV components of an image
FloatImage bilaYUV(const FloatImage &im, float sigmaRange, float sigmaY, float sigmaUV, float truncateDomain, bool clamp){

    // convert from RGB to YUV
    int size = 1 + 2 * ceil(max(sigmaY,sigmaUV) * truncateDomain);  // Use the larger sigma 
    float gaussian_d = 0;
    float gaussian_r = 0;
    float overall_val = 0;
    float k_normalize = 0 ;
    float diff = 0; 

    FloatImage rgb(im);
    FloatImage yuv(im.width(), im.height(), im.channels());
    FloatImage yuv_inter(im.width(), im.height(), im.channels());  // Make a copy, in case filtering the calculated pixels.
    yuv = rgb2yuv (rgb);
    yuv_inter = yuv; 
    // denoise Y and UV channels using different domain sigmas
    for (int i = 0; i < yuv_inter.width(); ++i)
        for (int j = 0; j < yuv_inter.height(); ++j)
            for (int k = 0; k < yuv_inter.channels(); ++k)
            {
                gaussian_d = 0;
                overall_val = 0;
                k_normalize = 0;
                for (int f = - size / 2; f <= size /2; ++f)
                    for (int m = - size /2; m <= size /2; ++m){
                        if (k == 0)
                            gaussian_d = exp(- ((f * f) + (m * m)) / (2 * sigmaY * sigmaY)) / (sigmaY * sigmaY * 2 * M_PI);
                        else
                            gaussian_d = exp(- ((f * f) + (m * m)) / (2 * sigmaUV * sigmaUV)) / (sigmaUV * sigmaUV * 2 * M_PI);

                        diff = yuv_inter.smartAccessor(i,j,k) - yuv_inter.smartAccessor(i + f,j + m,k);
                        gaussian_r = exp(- (diff * diff) / (2 * sigmaRange * sigmaRange)) / (sigmaRange * sqrt(2 * M_PI));
                        overall_val += gaussian_d * gaussian_r * yuv_inter.smartAccessor(i + f,j + m,k);  // No need for flipping, Gaussian is symmetric.
                        k_normalize += gaussian_d * gaussian_r;
                    }
                yuv(i,j,k) = overall_val / k_normalize;
            }

    // convert from YUV back to RGB
    rgb = yuv2rgb(yuv);
    
    return rgb;//  CHANGEME
}


/**************************************************************
 //                 FILTER CLASS FUNCTIONS                  //
 *************************************************************/


// write a convolution function for the filter class
FloatImage Filter::Convolve(const FloatImage &im, bool clamp) const {
    
    FloatImage imFilter(im.width(), im.height(), im.channels());
    float itvl;
    // implement convultion
    // Hint: use use Filter::operator()(x, y) to access (x,y) kernel location
   for (int i = 0; i < imFilter.width(); ++i)
        for (int j = 0; j < imFilter.height(); ++j)
            for (int p = 0; p < imFilter.channels(); ++p)
            {
                itvl = 0;
                for (int f = - this->width /2; f <= this->width /2; f++)
                    for (int m = - this->height /2; m <= this->height /2 ; m++){
                        itvl +=  im.smartAccessor(i - f,j - m,p,clamp) * Filter::operator()(this->width /2 + f,this->height /2 + m);      
                    }
                imFilter(i,j,p) = itvl;  
            }
    return imFilter;
}


FloatImage median_filter(const FloatImage &im){
    FloatImage output(im);
    FloatImage output_inter(im);
    std::vector<float> sort_val;


    for (int i = 0; i < im.width() -1; ++i)
        for (int j = 0; j < im.height() -1 ; ++j)
            for (int k = 0; k < im.channels(); ++k){
                sort_val.clear();
                // cout<<i<<"j = "<<j<<k<<endl;
                for (int m = -1; m <= 1; ++m)
                    for (int n = -1; n <= 1; ++n)
                        sort_val.push_back(output_inter.smartAccessor(i+m,j+n,k));
                std::sort(sort_val.begin(), sort_val.end());
                output(i,j,k) = sort_val[4];
                // cout<< output(i,j,k) <<endl;
            }
    return output;
}

/**************************************************************
 //               DON'T EDIT BELOW THIS LINE                //
 *************************************************************/

// Create an image of 0's with a value of 1 in the middle. This function
// can be used to test that you have properly set the kernel values in your
// Filter object. Make sure to set k to be larger than the size of your kernel
FloatImage impulseImg(const int &k){
    
    // initlize a kxkx1 image of all 0's
    FloatImage impulse(k, k, 1);
    
    // set the center pixel to have intensity 1
    int center = floor(k/2);
    impulse(center,center,0) = 1;
    
    return impulse;
}

Filter::Filter(const vector<float> &fData, const int &fWidth, const int &fHeight) {
    
    // TODO: check that width*heigh = length of filterVals and that width and height are odd
    if (fWidth*fHeight != fData.size() )
        throw OutOfBoundsException();
    if (fWidth %2 ==0 || fHeight %2 == 0)
        throw OutOfBoundsException();
    kernel = fData;
    width = fWidth;
    height = fHeight;

    
}

Filter::Filter(const int &fWidth, const int &fHeight) {
  width = fWidth;
  height = fHeight;
  kernel = std::vector<float>(width*height,0);
}

const float & Filter::operator()(int x, int y) const {
    if (x < 0 || x >= width)
        throw OutOfBoundsException();
    if ( y < 0 || y >= height)
        throw OutOfBoundsException();
    
    return kernel[x +y*width];
    
}


float & Filter::operator()(int x, int y) {
    if (x < 0 || x >= width)
        throw OutOfBoundsException();
    if ( y < 0 || y >= height)
        throw OutOfBoundsException();
    
    return kernel[x +y*width];
}
Filter::~Filter() {}
