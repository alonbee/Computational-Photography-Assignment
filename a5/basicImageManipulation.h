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

#ifndef __BASICIMAGEMANIPULATION__H
#define __BASICIMAGEMANIPULATION__H

#include "floatimage.h"

// Functions from PA02
FloatImage changeGamma(const FloatImage & im, float old_gamma, float new_gamma);
FloatImage exposure(const FloatImage &im, float factor);
static const float weight_init[3] = {0.299, 0.587, 0.114};
FloatImage color2gray(const FloatImage &im, const std::vector<float> &weights = std::vector<float>(weight_init, weight_init+3));
std::vector<FloatImage> lumiChromi(const FloatImage &im);
FloatImage brightnessContrastLumi(const FloatImage &im, float brightF, float contrastF, float midpoint = 0.3);
FloatImage rgb2yuv(const FloatImage &im);
FloatImage yuv2rgb(const FloatImage &im);
FloatImage saturate(const FloatImage &im, const float &k);
std::vector<FloatImage> spanish(const FloatImage &im);
FloatImage grayworld(const FloatImage & in);
FloatImage lumiChromi2rgb(const FloatImage &lumi, const FloatImage &chromi);

// Functions from PA00
FloatImage brightness(const FloatImage &im, float factor);
FloatImage contrast(const FloatImage &im, float factor, float midpoint = 0.5);


#endif
