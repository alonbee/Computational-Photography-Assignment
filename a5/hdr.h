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
#ifndef __hdr__h
#define __hdr__h

// hdr.h
// Assignment 5

#include "floatimage.h"
#include "basicImageManipulation.h"
#include <iostream>
#include <math.h>

// HDR
FloatImage computeWeight(const FloatImage &im, float epsilonMini=0.002, float epsilonMaxi=0.99);
float computeFactor(const FloatImage &im1, const FloatImage &w1, const FloatImage &im2, const FloatImage &w2);
FloatImage makeHDR(std::vector<FloatImage> &imSeq, float epsilonMini=0.002, float epsilonMaxi=0.99);

// Tone Mapping
FloatImage toneMap(const FloatImage &im, float targetBase=100, float detailAmp=3, bool useBila=false, float sigmaRange=0.1);
FloatImage exp10FloatImage(const FloatImage &im);
FloatImage log10FloatImage(const FloatImage &im);
float image_minnonzero(const FloatImage &im);

#endif
