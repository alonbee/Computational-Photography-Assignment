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
    // FloatImage im_chrominance(im.width(), im.height(), im.channels);
    // FloatImage im_luminance(im.width(), im.height(), 1);
    // std::vector<FloatImage> output;
    // Create luminance and chrominance images
    // output.push_back(im_luminance);
    // output.push_back(im_chrominance);
    // return output;
    FloatImage im_chrominance = im;
    FloatImage im_luminance(im.width(), im.height(), 1);
    std::vector<FloatImage> output;
    // Create luminance and chrominance images
    im_luminance = color2gray(im);
    for (int i = 0 ; i < im.width(); i++) {
        for (int j = 0 ; j < im.height(); j++) {
            for (int c = 0 ; c < im.channels(); c++ ){
                im_chrominance(i,j,c) = im_chrominance(i,j,c) / im_luminance(i,j);
            }
        }
    }
    output.push_back(im_luminance);
    output.push_back(im_chrominance);
    return output;
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
