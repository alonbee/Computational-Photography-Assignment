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

#include "basicImageManipulation.h"
#include <iostream>
#include <math.h>
        
using namespace std;

// create a new image that is k times bigger than the input by using nearest neighbor interpolation
FloatImage scaleNN(const FloatImage &im, float factor)
{
	// create new FloatImage that is factor times larger than input image
    FloatImage output(floor(im.width() * factor), floor(im.height() * factor), im.channels());
	// loop over new FloatImage pixels and set appropriate values (smartAccessor is probably overkill here)
    for (int i = 0; i < output.width();++i)
        for(int j = 0; j < output.height();++j)
            for (int k = 0; k < output.channels(); ++k)
            {
                output(i,j,k) = im(round(i / factor), round(j / factor) , k);
            }
    return output; // CHANGE ME
}



float interpolateLin(const FloatImage &im, float x, float y, int z, bool clamp)
{
    float final = 0;
    float L = 0;
    float U = 0;

    L = (1-(x - floor(x))) * im.smartAccessor(floor(x),floor(y),z,clamp) + (x - floor(x)) * im.smartAccessor(ceil(x),floor(y),z,clamp);
    U = (1-(x - floor(x))) * im.smartAccessor(floor(x),ceil(y),z,clamp) + (x - floor(x)) * im.smartAccessor(ceil(x),ceil(y),z,clamp);
    
    final = (1-(y - floor(y)))* L + (y - floor(y)) * U;
    
    // return final float value
    return final; // CHANGE ME
}


// create a new image that is k times bigger than the input by using bilinear interpolation
FloatImage scaleLin(const FloatImage &im, float factor) // Todo: should it have block effect?
{
	// create new FloatImage that is factor times larger than input image
    FloatImage output(floor(im.width() * factor), floor(im.height() * factor), im.channels());
    // loop over new FloatImage pixels and set appropriate values (use interpolateLin())
    if (factor < 1)
        throw parameteroutofbound();


    for (int i = 0; i < output.width(); ++i)
        for (int j = 0; j < output.height(); ++j)
            for (int k = 0; k < output.channels(); ++k)
        {
            output(i,j,k) = interpolateLin(im, float(i / factor), float (j / factor), k, true) ;
        }

    return output; // CHANGE ME
}


// rotate an image around its center by theta
FloatImage rotate(const FloatImage &im, float theta)  //TODO: Not perfect
{
    //rotate im around its center by theta
    // int diagonal = ceil(sqrt(im.width()^2 + im.height()^2));
    // FloatImage output(diagonal,diagonal,im.channels());
    FloatImage output(im.width(),im.height(),im.channels());
    std::vector<int> center_org(2,0);
    // std::vector<int> center_now(2,0);

    center_org[0] = im.width() / 2 -1;
    center_org[1] = im.height() / 2 -1;

    // center_now[0] = output.width() / 2;
    // center_now[1] = output.height() /2; 
    float len;
    float theta_now;
    float theta_orignal;
    float delta_x;
    float delta_y;
    float x;
    float y;

    for(int i = 0; i < output.size(); ++i)
        output(i) = 0;

    for (int i = 0; i < output.width(); ++i)
        for(int j = 0; j < output.height();++j)
            for(int k = 0; k < output.channels();++k){
                delta_x = i - center_org[0];
                delta_y = j - center_org[1];
                len = sqrt(delta_x * delta_x + delta_y * delta_y);

                if(delta_x != 0){
                    theta_now = atan(delta_y / delta_x);
                    if ( delta_x < 0)
                        theta_now += M_PI;
                }
                else if (delta_y > 0)
                    theta_now = M_PI / 2;
                else if (delta_y < 0)
                    theta_now = - M_PI / 2; 
                else if (delta_y == 0)
                    theta_now = 0;  

                theta_orignal = theta_now + theta;
                x = center_org[0] + len * cos(theta_orignal);
                y = center_org[1] + len * sin(theta_orignal);
                // if (j== center_org[1] || j == center_org[1] + 1 )
                // cout<< "i = "<<i<<"  j= "<<j<<" x = "<<x<<" y ="<<y<<" len = "<<len<<"theta_orignal"<<theta_orignal<<" "<<"theta_now"<<theta_now<<endl;

                 if (x<output.width() && x >=0 && y<output.height() && y >=0)  // Only manipulate the necessary pixels
                 output(i,j,k) = interpolateLin(im, x, y, k, true) ;  
            }
    return output; // CHANGE ME
}

/**************************************************************
 //               DON'T EDIT BELOW THIS LINE                //
 //          Functions from previous assignments            //
 *************************************************************/

// The implementation of these functions will be given to you after the pset 0 deadline
FloatImage brightness(const FloatImage &im, float factor) {
    return im * factor;
}

FloatImage contrast(const FloatImage &im, float factor, float midpoint) {
    return (im - midpoint)*factor+midpoint;
}



FloatImage changeGamma(const FloatImage & im, float old_gamma, float new_gamma) {
    // FloatImage output(im.width(), im.height(), im.channels());
    // Figure out what power to take the values of im, to get the values of output
    // return output;
    float exponent = new_gamma/old_gamma;
    FloatImage output = im;
    for (int i = 0 ; i < im.size();i++)
        output(i) = pow(im(i), exponent);
    return output;
}

// Change the exposure of the image. This is different than brightness because
// it means you multiply an image's intensity in the linear domain.
FloatImage exposure(const FloatImage &im, float factor) {
    FloatImage output = changeGamma(im, 1/2.2, 1.0);
    output = output*factor;
    output = changeGamma(output, 1.0, 1/2.2);
    return output;
    
}


FloatImage color2gray(const FloatImage &im, const std::vector<float> &weights) {
    // FloatImage output(im.width(), im.height(), 1);
    // Convert to grayscale
    // return output;
    float normalization = 0;
    for (unsigned int i = 0 ; i < weights.size();i++) {
      normalization += weights[i];
    }
    FloatImage output(im.width(), im.height(), 1);
    for (int i = 0 ; i < im.width(); i++ ) {
        for (int j = 0 ; j < im.height(); j++ ) {
            output(i,j,0) = im(i,j,0) * weights[0] + im(i,j,1) * weights[1] + im(i,j,2) *weights[2];
        }
    }
    output = output/normalization;
    return output;
}


// For this function, we want two outputs, a single channel luminance image
// and a three channel chrominance image. Return them in a vector with luminance first
std::vector<FloatImage> lumiChromi(const FloatImage &im) {
    // Create luminance and chrominance images
    FloatImage im_luminance = color2gray(im);
    FloatImage im_chrominance = im/im_luminance;
    std::vector<FloatImage> output;
    output.push_back(im_luminance);
    output.push_back(im_chrominance);
    return output;
}


// go from a luminance/crominance images back to an rgb image
FloatImage lumiChromi2rgb(const FloatImage &lumi, const FloatImage &chromi) {
    return chromi * lumi;
}

// Modify brightness then contrast
FloatImage brightnessContrastLumi(const FloatImage &im, float brightF, float contrastF, float midpoint) {
    // Create output image of appropriate size
    // Modify brightness, then contrast of luminance image
    // return output;
    std::vector<FloatImage> output = lumiChromi(im);
    FloatImage im_luminance = output[0];
    FloatImage im_chrominance = output[1];
    im_luminance = brightness(im_luminance, brightF);
    im_chrominance = contrast(im_chrominance, contrastF, midpoint);
    for (int i = 0 ; i < im.width(); i++ ){
        for (int j = 0 ; j < im.height(); j++) {
            for (int c = 0; c < im.channels(); c++) {
                im_chrominance(i,j,c) = im_chrominance(i,j,c) * im_luminance(i,j);
            }
        }
    }
    return im_chrominance;
}

FloatImage rgb2yuv(const FloatImage &im) {
    // Create output image of appropriate size
    // Change colorspae
    // return output;
    FloatImage output(im.width(), im.height(), im.channels());
    for (int i = 0; i < im.width(); i++) {
        for (int j = 0 ; j < im.height(); j++) {
            output(i,j,0) = 0.299 * im(i,j,0) + 0.587 * im(i,j,1) + 0.114 * im(i,j,2);
            output(i,j,1) = -0.147 * im(i,j,0) + (-0.289) *im(i,j,1) + 0.436 * im(i,j,2);
            output(i,j,2) = 0.615 * im(i,j,0) + -0.515 * im(i,j,1) + (-0.100) * im(i,j,2);
        }
    }
    return output;
}

FloatImage yuv2rgb(const FloatImage &im) {
    // Create output image of appropriate size
    // Change colorspae
    // return output;
    FloatImage output(im.width(), im.height(), im.channels());
    for (int i = 0; i < im.width(); i++) {
        for (int j = 0 ; j < im.height(); j++) {
            output(i,j,0) =  im(i,j,0) + 0 * im(i,j,1) + 1.14 * im(i,j,2);
            output(i,j,1) =  im(i,j,0) -0.395 *im(i,j,1) -0.581 * im(i,j,2);
            output(i,j,2) =  im(i,j,0) + 2.032 * im(i,j,1) + 0 * im(i,j,2);
        }
    }
    return output;
}

FloatImage saturate(const FloatImage &im, float factor) {
    // Create output image of appropriate size
    // Saturate image
    // return output;
    FloatImage output = rgb2yuv(im);
    for (int i = 0 ; i < im.width(); i++) {
        for (int j = 0 ; j < im.height(); j++) {
            output(i,j,1) = output(i,j,1) * factor;
            output(i,j,2) = output(i,j,2) * factor;
        }
    }
    output = yuv2rgb(output);
    return output;
}

// Return two images in a C++ vector
std::vector<FloatImage> spanish(const FloatImage &im) {
    // Remember to create the output images and the output vector
    // Push the images onto the vector
    // Do all the required processing
    // Return the vector, color image first
    FloatImage output_L = color2gray(im);
    FloatImage output_C = rgb2yuv(im);
    
    for (int i = 0; i < im.width(); i++) {
        for (int j = 0; j < im.height(); j++) {
            output_C(i,j,0) = 0.5;
            output_C(i,j,1) = -output_C(i,j,1);
            output_C(i,j,2) = -output_C(i,j,2);
        }
    }
    output_C = yuv2rgb(output_C);
    int bdot1 = floor(im.width()/2);
    int bdot2 = floor(im.height()/2);
    
    output_L(bdot1, bdot2,0) = 0;
    output_C(bdot1,bdot2,0) = 0;
    output_C(bdot1,bdot2,1) = 0;
    output_C(bdot1,bdot2,2) = 0;
    
    std::vector<FloatImage> output;
    output.push_back(output_C);
    output.push_back(output_L);
    return output;
}


// White balances an image using the gray world assumption
FloatImage grayworld(const FloatImage & im) {
    // Your code goes here
    float kr = 0, kg = 0, kb = 0;
    float N = im.width()*im.height();
    for (int i = 0 ; i < im.width(); i++) {
        for (int j = 0 ; j < im.height(); j++) {
            kr = kr + im(i,j,0);
            kg = kg + im(i,j,1);
            kb = kb + im(i,j,2);
        }
    }
    kr = kr/N;
    kg = kg/N;
    kb = kb/N;
    float redFactor = kg / kr;
    float blueFactor = kg / kb;
    FloatImage output = im;
    for (int i = 0 ; i < im.width(); i++) {
        for (int j = 0 ; j < im.height();j ++) {
      output(i,j,0) = output(i,j,0)*redFactor;
      output(i,j,1) = output(i,j,1);
      output(i,j,2) = output(i,j,2)*blueFactor;
        }
    }
    return output;
}

