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

#include "hdr.h"

using namespace std;


// test the compute weight function
void testComputeWeight() {

    // load in image and invert gamma correction
    FloatImage im1("./Input/ante1-1.png");
    im1 = changeGamma(im1, 1.0/2.2, 1.0f);
    
    // compute the weight image and save it
    FloatImage weight = computeWeight(im1, 0.01, 0.99);
    weight.write("./Output/weight.png");

}

// test the compute factor function
void testComputeFactor() {

    // load in 2 images
    FloatImage im1("./Input/ante2-1.png");
    FloatImage im2("./Input/ante2-2.png");
    
    // invert gamma correction
    im1 = changeGamma(im1, 1.0/2.2, 1.0f);
    im2 = changeGamma(im2, 1.0/2.2, 1.0f);

    // compute weight images and save them
    FloatImage w1 = computeWeight(im1);
    FloatImage w2 = computeWeight(im2);
    w1.write("./Output/ante2-w1.png");
    w2.write("./Output/ante2-w2.png");

    // compute the factor
    float factor = computeFactor(im1, w1, im2, w2);
    cout << "factor: " << factor << endl;

}

// test the hdr merging function
void testMakeHDR() {

    // load an image sequence
    vector<FloatImage> imSeq;
    imSeq.push_back(FloatImage("./Input/design-1.png"));
    imSeq.push_back(FloatImage("./Input/design-2.png"));
    imSeq.push_back(FloatImage("./Input/design-3.png"));
    imSeq.push_back(FloatImage("./Input/design-4.png"));
    imSeq.push_back(FloatImage("./Input/design-5.png"));
    imSeq.push_back(FloatImage("./Input/design-6.png"));
    imSeq.push_back(FloatImage("./Input/design-7.png"));

    // generate an hdr image
    FloatImage hdr = makeHDR(imSeq);

    // save out images clipped to different ranges.
    float maxVal = hdr.max();
    FloatImage hdrScale0 = hdr/maxVal;
    hdrScale0.write("./Output/scaledHDR_design_0.png");
    FloatImage hdrScale2 = (2e2)*hdr/maxVal;
    hdrScale2.write("./Output/scaledHDR_design_2.png");
    FloatImage hdrScale4 = (2e4)*hdr/maxVal;
    hdrScale4.write("./Output/scaledHDR_design_4.png");
    FloatImage hdrScale6 = (2e6)*hdr/maxVal;
    hdrScale6.write("./Output/scaledHDR_design_6.png");
    FloatImage hdrScale8 = (2e8)*hdr/maxVal;
    hdrScale8.write("./Output/scaledHDR_design_8.png");
    FloatImage hdrScale10 = (2e10)*hdr/maxVal;
    hdrScale10.write("./Output/scaledHDR_design_10.png");

}



// HDR and Tone Mapping on Ante2 images
void testToneMapping_ante2() {
    
    // load images
    vector<FloatImage> imSeq;
    imSeq.push_back(FloatImage("./Input/ante2-1.png"));
    imSeq.push_back(FloatImage("./Input/ante2-2.png"));
    
    // create hdr image
    FloatImage hdr = makeHDR(imSeq);
    
    // tone map with gaussian blur
    toneMap(hdr, 100, 1, false).write("./Output/ante2-tonedHDRsimple-gauss.png");
    // tone map with bilaterial
    toneMap(hdr, 100, 3, true, 0.1).write("./Output/ante2-tonedHDRsimple-bilateral.png");
    
}

// HDR and Tone Mapping on Ante3 images
void testToneMapping_ante3() {
    
    // load images
    vector<FloatImage> imSeq;
    imSeq.push_back(FloatImage("./Input/ante3-1.png"));
    imSeq.push_back(FloatImage("./Input/ante3-2.png"));
    imSeq.push_back(FloatImage("./Input/ante3-3.png"));
    imSeq.push_back(FloatImage("./Input/ante3-4.png"));

    // create hdr image
    FloatImage hdr = makeHDR(imSeq);

    // tone map with gaussian blur
    toneMap(hdr, 100, 1).write("./Output/ante3-tonedHDRsimple-gauss.png");
    // tone map with bilaterial
    toneMap(hdr, 100, 3, true, 0.1).write("./Output/ante3-tonedHDRsimple-bilateral.png");

}

// HDR and Tone Mapping on Design Images
void testToneMapping_design() {
    
    // load images
    vector<FloatImage> imSeq;
    imSeq.push_back(FloatImage("./Input/design-1.png"));
    imSeq.push_back(FloatImage("./Input/design-2.png"));
    imSeq.push_back(FloatImage("./Input/design-3.png"));
    imSeq.push_back(FloatImage("./Input/design-4.png"));
    imSeq.push_back(FloatImage("./Input/design-5.png"));
    imSeq.push_back(FloatImage("./Input/design-6.png"));
    imSeq.push_back(FloatImage("./Input/design-7.png"));
    
    // create hdr image
    FloatImage hdr = makeHDR(imSeq);
    
    // tone map with gaussian blur
    toneMap(hdr, 100, 1, false).write("./Output/design-tonedHDRsimple-gauss.png");
    
    // Note: bilaterial filtering these images takes a very long time. It is not
    // necessary to attempt this for testing
    // toneMap(hdr, 100, 3, true, 0.1).write("./Output/design-tonedHDRsimple-bilateral.png");
}

// HDR and Tone Mapping on Boston Images
void testToneMapping_boston() {
    
    // load images
    vector<FloatImage> imSeq;
    imSeq.push_back(FloatImage("./Input/boston-1.png"));
    imSeq.push_back(FloatImage("./Input/boston-2.png"));
    imSeq.push_back(FloatImage("./Input/boston-3.png"));
    
    // create hdr image
    FloatImage hdr = makeHDR(imSeq);
    
    // tone map with gaussian blur
    toneMap(hdr, 100, 1, false).write("./Output/boston-tonedHDRsimple-gauss.png");
    // tone map with bilaterial
    toneMap(hdr, 100, 3, true, 0.1).write("./Output/boston-tonedHDRsimple-bilateral.png");

}

void testToneMapping_other() {
    
    // load images
    vector<FloatImage> imSeq;
    imSeq.push_back(FloatImage("./Input/nyc-1.png"));
    imSeq.push_back(FloatImage("./Input/nyc-2.png"));
    
    // create hdr image
    FloatImage hdr = makeHDR(imSeq);
    
    // tone map with gaussian blur
    toneMap(hdr, 100, 1, false).write("./Output/nyc-tonedHDRsimple-gauss.png");
    // tone map with bilaterial
    toneMap(hdr, 100, 3, true, 0.1).write("./Output/nyc-tonedHDRsimple-bilateral.png");

}

void testmyownpic() {
    
    // load images
    vector<FloatImage> imSeq;
    imSeq.push_back(FloatImage("./Input/amine_1.png"));
    imSeq.push_back(FloatImage("./Input/amine_2.png"));
    imSeq.push_back(FloatImage("./Input/amine_3.png"));
    imSeq.push_back(FloatImage("./Input/amine_4.png"));
    imSeq.push_back(FloatImage("./Input/amine_5.png"));
    imSeq.push_back(FloatImage("./Input/amine_6.png"));
    
    // create hdr image
    FloatImage hdr = makeHDR(imSeq);
      cout<<"HDR finished" <<endl;

    // tone map with gaussian blur
    // toneMap(hdr, 100, 1, false).write("./Output/mine-tonedHDRsimple-gauss.png");
    // tone map with bilaterial
    toneMap(hdr, 100, 3, true, 0.1).write("./Output/mine-tonedHDRsimple-bilateral.png");

}
// This is a way for you to test your functions.
// We will not grade the contents of the main function
int main() {
    
    
    // uncomment these test functions as you complete the assignment
    
    // create an hdr image
    // testComputeWeight();
    // testComputeFactor();
    // testMakeHDR();
    
    // tone map an hdr image
    // testToneMapping_ante2();
    // testToneMapping_ante3();
    // testToneMapping_boston();
    // testToneMapping_design();
    testmyownpic();
        // testToneMapping_other();
}
