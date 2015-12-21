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

#ifndef __ALIGN__H
#define __ALIGN__H

#include "floatimage.h"
#include <iostream>
#include <math.h>

FloatImage denoiseSeq(const std::vector<FloatImage> &imgs);
FloatImage logSNR(const std::vector<FloatImage> &imSeq, float scale=1.0/20.0);
std::vector<int> align(const FloatImage &im1, const FloatImage &im2, int maxOffset=20);
FloatImage alignAndDenoise(const std::vector<FloatImage> &imSeq, int maxOffset=20);
FloatImage split(const FloatImage &sergeyImg);
FloatImage sergeyRGB(const FloatImage &sergeyImg, int maxOffset=20);

FloatImage roll(const FloatImage &im, int xRoll, int yRoll); 
 
#endif
