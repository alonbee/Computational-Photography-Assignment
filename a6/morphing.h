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

#ifndef __morphing__h
#define __morphing__h

// morphing.h
// Assignment 6

#include "floatimage.h"
#include "basicImageManipulation.h"
#include "exceptions.h"
#include <iostream>
#include <math.h>


    
class Segment
{
public:
    
    std::vector<float> P;
    std::vector<float> Q;
    std::vector<float> PQ;
    std::vector<float> perpPQ;
    std::vector<float> PQDivByPQlength2;
    std::vector<float> perpPQDivByPQlength;
    float PQ2;
    float PQlength;
    
    float getU(float x, float y);
    float getV(float x, float y);
    float dist(float x, float y);
    std::vector<float> UVtoX(float u, float v);
    float weight(float x, float y, float a, float b, float p);

    
    //Constructor
    Segment(float x1, float y1, float x2, float y2);
    
    // Destructor. Because there is no explicit memory management here, this doesn't do anything
    ~Segment();

    
// The following are static helper functions
    
    static float dot(const std::vector<float> &vec1, const std::vector<float> &vec2);
    static std::vector<float> subtract(const std::vector<float> &vec1, const std::vector<float> &vec2);
    static std::vector<float> scalarMult(const std::vector<float> &vec, float factor);

};

FloatImage warpBy1(const FloatImage &im, Segment &segBefore, Segment &segAfter); 
FloatImage warp(const FloatImage &im, std::vector<Segment> &segsBefore, std::vector<Segment> &segsAfter, float a=10.0, float b=1.0, float p=1.0);
std::vector<FloatImage> morph(const FloatImage &im1, const FloatImage &im2, std::vector<Segment> &segsBefore, std::vector<Segment> &segsAfter, int N=1, float a=10.0, float b=1.0, float p=1.0);

 
#endif
