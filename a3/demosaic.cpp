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

#include "demosaic.h"
#include "utils.h"
#include <math.h>

using namespace std;


// takes as input a raw image and returns a single-channel
// 2D image corresponding to the green channel using simple interpolation
FloatImage basicGreen(const FloatImage &raw, int offset)
{
	// create an appropriately sized output FloatImage
	FloatImage output(raw);	//	CHANGEME
	// cout << "num of channels" <<output.channels() <<endl;
	
	if(offset == 1){
		for (int x = 2; x < output.width() - 2; ++x)
			for (int y = 2; y < output.height() -2; ++y)
				if ((x%2 == 0 && y%2 ==0) || (x%2 == 1 && y%2 ==1))
				output(x,y,1) = (output(x-1,y,1) + output(x+1,y,1) + output(x,y+1,1) + output(x,y-1,1))/4;	
	}
		else if(offset == 0){
			for (int x = 2; x < output.width() - 2; ++x)
				for (int y = 2; y < output.height() -2; ++y)
					if ((x%2 == 1 && y%2 ==0) || (x%2 == 0 && y%2 ==1))
					output(x,y,1) = (output(x-1,y,1) + output(x+1,y,1) + output(x,y+1,1) + output(x,y-1,1))/4;	
		}
		else 
			throw OffsetException();
	// perform a simple interpolation given the offset
    return output;
}              

// takes as input a raw image and returns a single-channel
// 2D image corresponding to the red or blue channel using simple interpolation
FloatImage basicRorB(const FloatImage &raw, int offsetX, int offsetY)
{
	cout << "RorB(" << offsetX << ", " << offsetY << ")" << endl;
	// create an appropriately sized output FloatImage
	FloatImage output(raw); //	CHANGEME
	if(offsetX > 1 || offsetY >1)
		throw OffsetException();
	// horizontal interpolation
	// vertical interpolation	
	// diagonal interpolation

	if (offsetX == 0 && offsetY == 0){
		for (int x = 2; x < output.width() - 2; ++x)
			for (int y = 2; y < output.height() -2; ++y){
				if (x%2 ==1 && y%2 ==0)
				output(x,y,2) = (output(x-1,y,2) + output(x+1,y,2)) / 2;	
				else if (x%2 == 0 && y%2 ==1)
					output(x,y,2) = (output(x,y-1,2) + output(x,y+1,2)) /2;
					else if (x%2 == 1 && y%2 == 1)
						output(x,y,2) = (output(x-1,y+1,2) + output(x+1,y+1,2) + output(x-1,y-1,2) + output(x+1,y-1,2)) / 4;
			}

	}
	else if (offsetX == 1 && offsetY == 1){	
		for (int x = 2; x < output.width() - 2; ++x)
			for (int y = 2; y < output.height() -2; ++y){
				if (x%2 ==1 && y%2 ==0)
				output(x,y,0) = (output(x,y-1,0) + output(x,y+1,0)) / 2;	
				else if (x%2 == 0 && y%2 ==1)
					output(x,y,0) = (output(x-1,y,0) + output(x+1,y,0)) /2;
					else if (x%2 == 0 && y%2 == 0)
						output(x,y,0) = (output(x-1,y+1,0) + output(x+1,y+1,0) + output(x-1,y-1,0) + output(x+1,y-1,0)) / 4;
		}
	}	
    return output;
}

// takes as input a raw image and returns an rgb image
// using simple interpolation to demosaic each of the channels
FloatImage basicDemosaic(const FloatImage &raw, int offsetGreen, int offsetRedX, int offsetRedY, int offsetBlueX, int offsetBlueY)
{	

	// create an appropriately sized output FloatImage
    FloatImage output(raw);	// CHANGEME

    //use basicGreen() and basicRorB() to demosaic the raw image
    output = basicGreen(output,offsetGreen);
    output = basicRorB(output,offsetRedX,offsetRedY);
    output = basicRorB(output,offsetBlueX,offsetBlueY);

	return output;
}

// takes a raw image and outputs a single-channel
// image corresponding to the green channel taking into account edges
FloatImage edgeBasedGreen(const FloatImage &raw, int offset)
{
	// create an appropriately sized output FloatImage
	FloatImage output(raw);	//	CHANGEME
	float hor_diff; 
	float ver_diff; 

	//	interpolate horizontally or vertically
	if(offset == 1){
		for (int x = 2; x < output.width() - 2; ++x)
			for (int y = 2; y < output.height() -2; ++y)
				if ((x%2 == 0 && y%2 ==0) || (x%2 == 1 && y%2 ==1)){
					hor_diff = fabs(output(x+1,y,1) - output(x-1,y,1));
					ver_diff = fabs(output(x,y+1,1) - output(x,y-1,1));
					if (ver_diff <= hor_diff)
						output(x,y,1) = (output(x,y+1,1) + output(x,y-1,1)) / 2;
						else
						output(x,y,1) = (output(x+1,y,1) + output(x-1,y,1)) / 2; 
				}			
	}	
		else if(offset == 0){
			for (int x = 2; x < output.width() - 2; ++x)
				for (int y = 2; y < output.height() -2; ++y)
					if ((x%2 == 1 && y%2 ==0) || (x%2 == 0 && y%2 ==1)){
					hor_diff = fabs(output(x+1,y,1) - output(x-1,y,1));
					ver_diff = fabs(output(x,y+1,1) - output(x,y-1,1));
					if (ver_diff <= hor_diff)
						output(x,y,1) = (output(x,y+1,1) + output(x,y-1,1)) / 2;
						else
						output(x,y,1) = (output(x+1,y,1) + output(x-1,y,1)) / 2;
					} 
		}
		else 
			throw OffsetException();

    return output;
}

// takes as input a raw image and returns an rgb image
// using edge-based green demosaicing for the green channel and
// simple interpolation to demosaic the red and blue channels
FloatImage edgeBasedGreenDemosaic(const FloatImage &raw, int offsetGreen, int offsetRedX, int offsetRedY, int offsetBlueX, int offsetBlueY)
{
	// create an appropriately sized output FloatImage
    FloatImage output(raw);	// CHANGEME

    //use edgeBasedGreen() and basicRorB() to demosaic the raw image
    output = edgeBasedGreen(output,offsetGreen);
    output = basicRorB(output,offsetRedX,offsetRedY);
    output = basicRorB(output,offsetBlueX,offsetBlueY);
	return output;
}


// takes as input a raw image and returns a single-channel
// 2D image corresponding to the red or blue channel using green based interpolation
FloatImage greenBasedRorB(const FloatImage &raw, FloatImage &green, int offsetX, int offsetY)
{
	// create an appropriately sized output FloatImage
    FloatImage output(raw);	// CHANGEME
    FloatImage green_base(green); 

    // like basicRorB(), but interpolate the differences R-G or B-G 
    if (offsetX == 0 && offsetY == 0){
		for (int x = 2; x < output.width() - 2; ++x)
			for (int y = 2; y < output.height() -2; ++y){
				//if ( x%2 == 0 && y%2 == 0 )  // ??
				output(x,y,2) = output(x,y,2) - green_base(x,y,1); }

				output = basicRorB(output,offsetX,offsetY);

		for (int x = 2; x < output.width() - 2; ++x)
			for (int y = 2; y < output.height() -2; ++y)
				output(x,y,2) = output(x,y,2) + green_base(x,y,1);

	}
	else if (offsetX == 1 && offsetY == 1){	
		for (int x = 2; x < output.width() - 2; ++x)
			for (int y = 2; y < output.height() -2; ++y){
				//if (x%2 == 1 && y%2 ==1)	//??
					output(x,y,0) = output(x,y,0) - green_base(x,y,1);
			}

			output = basicRorB(output,offsetX,offsetY);

		for (int x = 2; x < output.width() - 2; ++x)
			for (int y = 2; y < output.height() -2; ++y)

				output(x,y,0) = output(x,y,0) + green_base(x,y,1);

		}
	else 
		throw OffsetException();

    return output;
}

// takes as input a raw image and returns an rgb image
// using edge-based green demosaicing for the green channel and
// simple green based demosaicing of the red and blue channels
FloatImage improvedDemosaic(const FloatImage &raw, int offsetGreen, int offsetRedX, int offsetRedY, int offsetBlueX, int offsetBlueY)
{
	// create an appropriately sized output FloatImage
    FloatImage output(raw);	// CHANGEME
	FloatImage green_base(raw);
    //use edgeBasedGreen() and greenBasedRorB() to demosaic the image

	green_base = edgeBasedGreen(green_base,offsetGreen);
	for (int x = 2; x < output.width() - 2; ++x)
		for (int y = 2; y < output.height() -2; ++y)
			output(x,y,1) = green_base(x,y,1);
	output = greenBasedRorB(output,green_base,0,0);
	output = greenBasedRorB(output,green_base,1,1);
	return output;
}
