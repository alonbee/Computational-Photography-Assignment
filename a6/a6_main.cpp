#include "morphing.h"


using namespace std;


// test the smart accessor function
void testSmartAccessor(){
    
    // load an image and create 2 images that will test the smart accessor
    const FloatImage input("./Input/bear.png");
    input.write("./Output/bear.png");

    FloatImage clampTrue(input.width(), input.height(), input.channels());
    FloatImage clampFalse(input.width(), input.height(), input.channels());
    
    for (int z=0; z<input.channels(); z++){
        for (int x=0; x<input.width(); x++){
            for (int y=0; y<input.height(); y++){
                
                // replace non-valid pixel values with the value of the nearest pixel
                clampTrue(x,y,z) = input.smartAccessor(x-10, y-10, z, true);
                // replace non-valid pixel values with black (value=0)
                clampFalse(x,y,z) = input.smartAccessor(x-10, y-10, z);
            }
        }
    }
    
    clampTrue.write("./Output/smartAccessor_clampTrue.png");
    clampFalse.write("./Output/smartAccessor_clampFalse.png");
    
}

// a function to test scaling
void testScaling(){
    
    // load in the image and print out original size
    float fact = 2.253;
    const FloatImage bostonim("./Input/BostonRainbow-crop-400.png");
    printf("Boston image is %dx%dx%d\n", bostonim.width(), bostonim.height(), bostonim.channels());
    
    // scale using NN interpolation and print the size of the new image
    FloatImage scaledNN = scaleNN(bostonim, fact);
    scaledNN.write("./Output/BostonRainbow-scaled-NN.png");
    printf("Scaled-NN image is %dx%dx%d\n", scaledNN.width(), scaledNN.height(), scaledNN.channels());
    
    // scale using bilinear interpolation and print the size of the new image
    FloatImage scaledLin = scaleLin(bostonim, fact);
    scaledLin.write("./Output/BostonRainbow-scaled-Lin.png");
    printf("Scaled-Lin image is %dx%dx%d\n", scaledLin.width(), scaledLin.height(), scaledLin.channels());
    
}

// a function to test rotation for those who have done it
void testRotation(){
    
    float theta = 3.141/4;
    
    const FloatImage bostonim("./Input/BostonRainbow-crop-400.png");
    
    FloatImage rot = rotate(bostonim, theta);
    rot.write("./Output/BostonRainbow-rotated.png");
    printf("Scaled-Lin image is %dx%dx%d\n", rot.width(), rot.height(), rot.channels());
    
}


// test warp by 1
void testWarpBy1(){
    
    FloatImage bearim("./Input/bear.png");
    FloatImage fredo("./Input/fredo.png");

    // define before and after segments
    Segment segBefore(0,0, 10,0);
    Segment segAfter(10, 10, 30, 15);

    
    FloatImage warp1 = warpBy1(bearim, segBefore, segAfter);
        FloatImage warp2 = warpBy1(fredo, segBefore, segAfter);

    warp1.write("./Output/warpBy1.png");
    warp2.write("./Output/warpBy2.png");

}


// a function to test your morphing function with the fredo and werewolf images
void testMorph(){
    
    // load the images
    FloatImage fredo("./Input/fredo.png");
    FloatImage werewolf("./Input/werewolf.png");
    
    // print out the size of the two images
    printf("Fredo image is %dx%dx%d\n", fredo.width(), fredo.height(), fredo.channels());
    printf("Werewolf image is %dx%dx%d\n", werewolf.width(), werewolf.height(), werewolf.channels());


    // define the segments corresponding to fredo (segsBefore) and the werewolf (segsAfter)
    vector<Segment> segsBefore, segsAfter;
    
    segsBefore.push_back(Segment(87, 128, 109, 127));
    segsBefore.push_back(Segment(143, 127, 162, 131));
    segsBefore.push_back(Segment(96, 197, 129, 190));
    segsBefore.push_back(Segment(118, 221, 132, 200));
    segsBefore.push_back(Segment(140, 238, 165, 170));
    segsBefore.push_back(Segment(71, 242, 44, 196));
    segsBefore.push_back(Segment(9, 46, 34, 14));
    segsBefore.push_back(Segment(175, 66, 154, 27));
    
    segsAfter.push_back(Segment(83, 114, 107, 109));
    segsAfter.push_back(Segment(139, 103, 157, 101));
    segsAfter.push_back(Segment(100, 170, 132, 151));
    segsAfter.push_back(Segment(125, 198, 145, 171));
    segsAfter.push_back(Segment(142, 196, 158, 139));
    segsAfter.push_back(Segment(96, 211, 75, 180));
    segsAfter.push_back(Segment(11, 41, 33, 7));
    segsAfter.push_back(Segment(175, 56, 155, 13));
    
    
    // create an image morphing between fredo and werewolf at time t=0.5
    vector<FloatImage> imMorph = morph(fredo, werewolf, segsBefore, segsAfter);
    cout<<"okkk"<<endl;
    // write out images
    char buffer [50];
    for ( int n=0; n< 3; n++){
        
        FloatImage im = imMorph[n];
        sprintf (buffer, "./Output/fredo_werewolf_morph_%d.png", n);
        
        im.write(buffer);
	}
}

void testwp(){
    
    // load the images
    FloatImage fredo("./Input/fredo.png");
    FloatImage werewolf("./Input/werewolf.png");
    
    // print out the size of the two images
    printf("Fredo image is %dx%dx%d\n", fredo.width(), fredo.height(), fredo.channels());
    printf("Werewolf image is %dx%dx%d\n", werewolf.width(), werewolf.height(), werewolf.channels());


    // define the segments corresponding to fredo (segsBefore) and the werewolf (segsAfter)
    vector<Segment> segsBefore, segsAfter;
    
    segsBefore.push_back(Segment(87, 128, 109, 127));
    segsBefore.push_back(Segment(143, 127, 162, 131));
    segsBefore.push_back(Segment(96, 197, 129, 190));
    segsBefore.push_back(Segment(118, 221, 132, 200));
    segsBefore.push_back(Segment(140, 238, 165, 170));
    segsBefore.push_back(Segment(71, 242, 44, 196));
    segsBefore.push_back(Segment(9, 46, 34, 14));
    segsBefore.push_back(Segment(175, 66, 154, 27));
    
    segsAfter.push_back(Segment(83, 114, 107, 109));
    segsAfter.push_back(Segment(139, 103, 157, 101));
    segsAfter.push_back(Segment(100, 170, 132, 151));
    segsAfter.push_back(Segment(125, 198, 145, 171));
    segsAfter.push_back(Segment(142, 196, 158, 139));
    segsAfter.push_back(Segment(96, 211, 75, 180));
    segsAfter.push_back(Segment(11, 41, 33, 7));
    segsAfter.push_back(Segment(175, 56, 155, 13));
    
    
    // create an image morphing between fredo and werewolf at time t=0.5
    FloatImage imMorph(fredo.width(),fredo.height(),fredo.channels());
     imMorph = warp(fredo, segsBefore, segsAfter);
     imMorph.write("./Output/warptest.png");
}

void classMorph() {
    FloatImage me("./Input/wcmm.png");
    FloatImage classmate("./Input/fyyr.png");
    
// vector<Segment> segsBefore;
// segsBefore.push_back(Segment(170, 223, 197, 231)); 
// segsBefore.push_back(Segment(200, 231, 235, 222)); 
// segsBefore.push_back(Segment(199, 161, 199, 203)); 
// segsBefore.push_back(Segment(144, 234, 174, 271)); 
// segsBefore.push_back(Segment(119, 157, 114, 99)); 
// segsBefore.push_back(Segment(125, 77, 200, 44)); 
// segsBefore.push_back(Segment(264, 56, 300, 138)); 
// segsBefore.push_back(Segment(268, 212, 249, 259)); 
// segsBefore.push_back(Segment(147, 161, 182, 161)); 
// segsBefore.push_back(Segment(217, 161, 255, 163)); 
// segsBefore.push_back(Segment(149, 140, 186, 145)); 
// segsBefore.push_back(Segment(211, 144, 258, 142)); 
// segsBefore.push_back(Segment(181, 277, 224, 280)); 

// vector<Segment> segsAfter;
// segsAfter.push_back(Segment(176, 239, 198, 239)); 
// segsAfter.push_back(Segment(203, 237, 226, 237)); 
// segsAfter.push_back(Segment(203, 161, 203, 203)); 
// segsAfter.push_back(Segment(136, 245, 172, 282)); 
// segsAfter.push_back(Segment(109, 162, 115, 90)); 
// segsAfter.push_back(Segment(125, 68, 208, 40)); 
// segsAfter.push_back(Segment(235, 47, 287, 149)); 
// segsAfter.push_back(Segment(267, 223, 247, 267)); 
// segsAfter.push_back(Segment(145, 170, 182, 168)); 
// segsAfter.push_back(Segment(216, 166, 259, 168)); 
// segsAfter.push_back(Segment(146, 151, 188, 148)); 
// segsAfter.push_back(Segment(213, 151, 257, 146)); 
// segsAfter.push_back(Segment(187, 289, 218, 289)); 

// vector<Segment> segsBefore;
// segsBefore.push_back(Segment(179, 256, 219, 259)); 
// segsBefore.push_back(Segment(190, 286, 106, 336)); 
// segsBefore.push_back(Segment(87, 345, 68, 386)); 
// segsBefore.push_back(Segment(69, 389, 117, 431)); 
// segsBefore.push_back(Segment(190, 164, 180, 219)); 
// segsBefore.push_back(Segment(149, 172, 179, 175)); 
// segsBefore.push_back(Segment(210, 171, 263, 165)); 
// segsBefore.push_back(Segment(149, 145, 175, 143)); 
// segsBefore.push_back(Segment(202, 138, 252, 138)); 
// segsBefore.push_back(Segment(225, 295, 286, 267)); 
// segsBefore.push_back(Segment(169, 247, 151, 202)); 
// segsBefore.push_back(Segment(275, 372, 362, 329)); 
// segsBefore.push_back(Segment(229, 315, 226, 353)); 

// vector<Segment> segsAfter;
// segsAfter.push_back(Segment(212, 211, 265, 215)); 
// segsAfter.push_back(Segment(210, 244, 111, 252)); 
// segsAfter.push_back(Segment(97, 263, 45, 332)); 
// segsAfter.push_back(Segment(42, 341, 102, 393)); 
// segsAfter.push_back(Segment(244, 108, 231, 166)); 
// segsAfter.push_back(Segment(180, 112, 222, 115)); 
// segsAfter.push_back(Segment(262, 117, 329, 116)); 
// segsAfter.push_back(Segment(187, 81, 227, 85)); 
// segsAfter.push_back(Segment(263, 80, 335, 79)); 
// segsAfter.push_back(Segment(256, 258, 314, 228)); 
// segsAfter.push_back(Segment(189, 197, 174, 152)); 
// segsAfter.push_back(Segment(275, 326, 336, 255)); 
// segsAfter.push_back(Segment(213, 258, 207, 321));  

vector<Segment> segsBefore;
segsBefore.push_back(Segment(120, 161, 163, 157)); 
segsBefore.push_back(Segment(204, 157, 241, 154)); 
segsBefore.push_back(Segment(187, 141, 199, 198)); 
segsBefore.push_back(Segment(163, 242, 219, 241)); 
segsBefore.push_back(Segment(147, 295, 95, 258)); 

vector<Segment> segsAfter;
segsAfter.push_back(Segment(115, 164, 164, 171)); 
segsAfter.push_back(Segment(194, 173, 227, 172)); 
segsAfter.push_back(Segment(186, 157, 192, 216)); 
segsAfter.push_back(Segment(156, 258, 192, 256)); 
segsAfter.push_back(Segment(135, 291, 74, 253)); 
 


    vector<FloatImage> imMorph = morph(me, classmate, segsBefore, segsAfter,40);
    // cout<<"okkk"<<endl;
    // write out images
    char buffer [imMorph.size()];
    for ( int n=0; n< imMorph.size(); n++){
        
        FloatImage im = imMorph[n];
        sprintf (buffer, "./Output/fyy3/fyy_all%d.png", n);
        
        im.write(buffer);
    }



}

// This is a way for you to test your functions.
// We will not grade the contents of the main function
int main() {
    // uncomment these test functions as you complete the assignment to test your code
    // try { testSmartAccessor();}   catch(...) {cout << "testSmartAccessor Failed!" << endl;}
    // try { testScaling(); }        catch(...) {cout << "testScaling Failed!" << endl;}
    // try { testRotation(); }       catch(...) {cout << "testRotation Failed!" << endl;}
    // try { testWarpBy1(); }        catch(...) {cout << "testWarpBy1 Failed!" << endl;}
    // try { testMorph(); }          catch(...) {cout << "testMorph Failed!" << endl;}
    classMorph();
}
