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

#ifndef __DEMOSAIC__H
#define __DEMOSAIC__H

#include "floatimage.h"
#include <iostream>
#include <math.h>

FloatImage basicGreen(const FloatImage &raw, int offset=1);
FloatImage basicRorB(const FloatImage &raw, int offsetX, int offsetY);
FloatImage basicDemosaic(const FloatImage &raw, int offsetGreen=1, int offsetRedX=1, int offsetRedY=1, int offsetBlueX=0, int offsetBlueY=0);
FloatImage edgeBasedGreen(const FloatImage &raw, int offset=1);
FloatImage edgeBasedGreenDemosaic(const FloatImage &raw, int offsetGreen=1, int offsetRedX=1, int offsetRedY=1, int offsetBlueX=0, int offsetBlueY=0);
FloatImage greenBasedRorB(const FloatImage &raw, FloatImage &green, int offsetX, int offsetY);
FloatImage improvedDemosaic(const FloatImage &raw, int offsetGreen=1, int offsetRedX=1, int offsetRedY=1, int offsetBlueX=0, int offsetBlueY=0);

 
#endif
