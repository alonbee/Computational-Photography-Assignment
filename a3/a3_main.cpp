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

#include "align.h"
#include "demosaic.h"
#include <ctime>

using namespace std;


// test denoising function by averaging function and visualizing SNR
void testDenoiseSeq()
{    
    // load an image sequence
    vector<FloatImage> imSeq;
    int nImages = 50;
    char baseName[] = "./Input/aligned-ISO3200/1D2N-iso3200-";
    char fileName [100];
    
    // add each image to the vector of images
    for (int i = 1; i <= nImages; i++)
    {
        sprintf(fileName, "%s%d.png", baseName, i);
        imSeq.push_back(FloatImage(fileName));
    }
    
    
    // denoise by averaging across all of the images
    FloatImage imDenoise = denoiseSeq(imSeq);
    imDenoise.write("./Output/denoisedImg.png");
    
    // estimate the SNR of the image sequence
    FloatImage SNRIm = logSNR(imSeq);
    SNRIm.write("./Output/snrImg.png"); 
}

// test your function to find the optimal offset between 2 images
void testOffset()
{    
    // load the images that are not aligned
    FloatImage im1("./Input/green/noise-small-1.png");
    FloatImage im2("./Input/green/noise-small-2.png");
    
    // get the best offset for these 2 images
    int maxOffset = 10;
    vector<int> offset = align(im1, im2, maxOffset);
    
    // print out the offset values that best aligns im2 to im1
    cout << "x offset: " << offset[0] << " y offset: " << offset[1] << endl;
}

// test your denoising function by aligned averaging
void testDenoiseShiftSeq()
{    
    // load an image sequence
    vector<FloatImage> imSeq;
    int nImages = 9;
    char baseName[] = "./Input/green/noise-small-";
    char fileName [100];
    
    // add each image to the vector of images
    for (int i = 1; i <= nImages; i++)
    {
        sprintf(fileName, "%s%d.png", baseName, i);
        imSeq.push_back(FloatImage(fileName));
    }
    
    // denoise using just pure averaging - this should result in a blurry output
    FloatImage imDenoiseBlur = denoiseSeq(imSeq);
    imDenoiseBlur.write("./Output/denoisedImg_Avg.png");
    
    // denoise by averaging aligned images
    FloatImage imDenoiseShift = alignAndDenoise(imSeq, 5);
    imDenoiseShift.write("./Output/denoisedImg_AlignedAvg.png");
    
}

// test your demosaicing functions that use simple interpolation
void testBasicDemosaic()
{
    // load raw image (notice that although this is a 3 channel image each channel
    // has the same value so you can think of it as having only 1 channel)
    FloatImage raw("./Input/raw/signs-small.png");
    
    // demosaic the green channel
    FloatImage green = basicGreen(raw);
    green.write("./Output/green_basic.png");
    
    // demosaic the red channel
    FloatImage red = basicRorB(raw, 1, 1);
    red.write("./Output/red_basic.png");
    
    // demosaic the blue channel
    FloatImage blue = basicRorB(raw, 0, 0);
    blue.write("./Output/blue_basic.png");
    
    // generate a full rgb image by demaosicing each of the channels
    // and concatinating them together
    FloatImage rgb = basicDemosaic(raw);
    rgb.write("./Output/rgb_basic.png");
}


// test your demosaicing functions that use the green channel to incoperate information
// about edges in your image
void testGreenEdgeDemosaic()
{
    // load raw image (notice that although this is a 3 channel image each channel
    // has the same value so you can think of it as having only 1 channel)
    FloatImage raw("./Input/raw/signs-small.png");
    
    // demosaic the green channel using the edge-based method
    FloatImage green = edgeBasedGreen(raw);
    green.write("./Output/green_edge.png");
    
    // generate the rgb image using edge based demosaicing algorithm for the
    // green channel but the naive red and blue demosaicing algorithm
    FloatImage rgb_greenEdge = edgeBasedGreenDemosaic(raw);
    rgb_greenEdge.write("./Output/rgb_greenEdge.png");
    
    // generate the rgb image using the edge based demosaicing algorithm
    // for the green, red, and blue channels
    FloatImage rgb_fullEdge = improvedDemosaic(raw);
    rgb_fullEdge.write("./Output/rgb_fullEdge.png");
    
}


// test your functions to generate an rgb image from the Sergey images
void testSergey()
{
    // load a sergey image
    FloatImage sergeyImg("./Input/Sergey/01880v_third.png");
    
    // split the grayscale image to generate a rgb image
    FloatImage rgb = split(sergeyImg);
    rgb.write("./Output/sergeyRGB.png");
    
    // split the grayscale image and align the channels to
    // generate a clean RGB image
    FloatImage rgbAlign = sergeyRGB(sergeyImg);
    rgbAlign.write("./Output/sergeyRGB_align.png");
}



// This is a way for you to test your functions.
// We will not grade the contents of the main function
int main()
{    
    // uncomment these test functions as you complete the assignment
     testDenoiseSeq();
     // testOffset();      // Segmentation fault: 11 ?? Why?
     // testDenoiseShiftSeq();
       // testBasicDemosaic();
       // testGreenEdgeDemosaic();
    // testSergey();
}
