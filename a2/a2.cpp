              
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

#include "a2.h"
#include <math.h>
#include "utils.h"
#include <assert.h>
#include <iostream>

using namespace std;

// The implementation of these functions will be given to you after the assignment 1 deadline
FloatImage brightness(const FloatImage &im, float factor)
{
    return im * factor;
}


FloatImage contrast(const FloatImage &im, float factor, float midpoint)
{
    return (im - midpoint)*factor+midpoint;
}


FloatImage changeGamma(const FloatImage & im, float old_gamma, float new_gamma)
{
    // create an appropriately sized output FloatImage
    FloatImage output(im);
    float ratio = new_gamma / old_gamma;

    // Figure out what power to take the values of im, to get the values of output
   
     for (int i = 0; i < output.size(); ++i)
        output(i) = pow(output(i),ratio);

    return output;
}


// Change the exposure of the image. This is different than brightness because
// it means you multiply an image's intensity in the linear domain.
FloatImage exposure(const FloatImage &im, float factor)
{
    // create an appropriately sized output FloatImage
    FloatImage output(im);
    float odgma =  1/2.2;
    output = changeGamma(output,odgma,1);
    output *= factor;
    output = changeGamma(output,1,odgma);
    return output;
}


FloatImage color2gray(const FloatImage &im, const vector<float> &weights)  //output should be 3 dimensional, the third is 1 
{
        // Convert to grayscale
    if (weights.size() != im.channels())    // random number: 0.837375 0.162625 6.64281e-10
        throw MismatchedDimensionsException();

    FloatImage output(im);  // CHANGEME
    float wt1 = weights[0] / (weights[0] + weights[1] + weights[2]);
    float wt2 = weights[1] / (weights[0] + weights[1] + weights[2]);
    float wt3 = weights[2] / (weights[0] + weights[1] + weights[2]);

    unsigned int wd = im.width();
    unsigned int ht = im.height();
    //unsigned int ch = im.channels();
    
    for (int x=0; x < wd; ++x)
        for (int y=0; y<ht; ++y){
            output(x,y,0) = output(x,y,0) * wt1 + output(x,y,1) * wt2 + output(x,y,2) * wt3;   // TODO LOW EFFICIENCY ?
            // output(x,y,1) = output(x,y,0);
            // output(x,y,2) = output(x,y,0);
        }

    output.resize(wd,ht,1);   // TODO: need to be fixed
    return output;
}


// For this function, we want two outputs, a single channel luminance image 
// and a three channel chrominance image. Return them in a vector with luminance first
vector<FloatImage> lumiChromi(const FloatImage &im)
{
    vector<FloatImage> output;
    FloatImage ele_1(im);
    FloatImage ele_2(im);
    ele_1 = color2gray(ele_1);  // Lumi of the first 

   // ele_2 /= ele_1;


    for (int i = 0; i < ele_2.width(); ++i)
        for (int j = 0; j < ele_2.height(); ++j)
            for (int k = 0; k < ele_2.channels(); ++k)
        {
            ele_2(i,j,k) = ele_2(i,j,k) / ele_1(i,j,0);     // Must have the same channels.
        }
    
    output.push_back(ele_1);
    output.push_back(ele_2);
    // Create luminance and chrominance images
    // push the luminance and chrominance images onto the vector

    return output;
}

// Modify brightness then contrast
FloatImage brightnessContrastLumi(const FloatImage &im, float brightF, float contrastF, float midpoint)
{
    // Create output image of appropriate size
    // Modify brightness, then contrast of luminance image
    FloatImage output(im);
    vector<FloatImage> dcps = lumiChromi(output);
    FloatImage lum(dcps[0]);
    FloatImage cts(dcps[1]);
    lum = contrast(lum,contrastF,midpoint);
    lum = brightness(lum,brightF);  // Not important:  brightness first or contrast first
    output = cts * lum;     // Must use cts * lum becuase lum has only one channel. 
    return output;
}

FloatImage rgb2yuv(const FloatImage &im)
{
    // Create output image of appropriate size
    FloatImage rgb(im);
    FloatImage yuv(im);
    for (int i = 0; i < yuv.width(); ++i)
        for (int j = 0; j < yuv.height(); ++j)
        {
            yuv(i,j,0) = rgb(i,j,0) * 0.299 + rgb(i,j,1) * 0.587 + rgb(i,j,2) * 0.114;
            yuv(i,j,1) = rgb(i,j,0) * (-0.147) + rgb(i,j,1) * (-0.289) + rgb(i,j,2) * 0.436;
            yuv(i,j,2) = rgb(i,j,0) * 0.615 + rgb(i,j,1) * (-0.515) + rgb(i,j,2) * (-0.1);
        }
        // Change colorspace
    return yuv; 
}

FloatImage yuv2rgb(const FloatImage &im)
{
    // Create output image of appropriate size
    // Change colorspace
    FloatImage yuv(im);
    FloatImage rgb(im);
    for (int i = 0; i < rgb.width(); ++i)
        for (int j = 0; j < rgb.height(); ++j)
        {
            rgb(i,j,0) = yuv(i,j,0) + yuv(i,j,1) * 0 + yuv(i,j,2) * 1.14;
            rgb(i,j,1) = yuv(i,j,0) + yuv(i,j,1) * (-0.395) + yuv(i,j,2) * (-0.581);
            rgb(i,j,2) = yuv(i,j,0) + yuv(i,j,1) * (2.032) + yuv(i,j,2) * 0;
        }
    return rgb; 
}

FloatImage saturate(const FloatImage &im, float factor)
{
    // Create output image of appropriate size
    FloatImage output(im);
    output = rgb2yuv(output); 
    for (int i = 0; i < output.width(); ++i)
        for (int j = 0; j < output.height(); ++j)
        {
            output(i,j,1) = output(i,j,1) * factor;
            output(i,j,2) = output(i,j,2) * factor;
        }  
    output = yuv2rgb(output);

    // saturate image

    return output; 
}

// Return two images in a C++ vector
vector<FloatImage> spanish(const FloatImage &im)
{
    // Remember to create the output images and the output vector
    vector<FloatImage> output;
    FloatImage gray(im);
    FloatImage yuv(im);
    unsigned int x  = floor(gray.width() / 2);
    unsigned int y = floor(gray.height() / 2);
    // Do all the required processing
    gray = rgb2yuv(gray);
    // Push the images onto the vector
    yuv = rgb2yuv(yuv);
    for (int i = 0; i < yuv.width(); ++i)
        for (int j = 0; j < yuv.height(); ++j)
        {
            yuv(i,j,0) = 0.5;
            yuv(i,j,1) = - yuv(i,j,1);
            yuv(i,j,2) = - yuv(i,j,2);
            gray(i,j,1) = 0;
            gray(i,j,2) = 0;
        }  
        gray = yuv2rgb(gray);
        yuv = yuv2rgb(yuv);
        gray(x,y,0) =  gray(x,y,1) = gray(x,y,2) = 0;
        gray(x+1,y,0) =  gray(x+1,y,1) = gray(x+1,y,2) = 0;
        gray(x,y+1,0) =  gray(x,y+1,1) = gray(x,y+1,2) = 0;
        gray(x+1,y+1,0) =  gray(x+1,y+1,1) = gray(x+1,y+1,2) = 0;

        yuv(x,y,0) =  yuv(x,y,1) = yuv(x,y,2) = 0;
        yuv(x+1,y,0) =  yuv(x+1,y,1) = yuv(x+1,y,2) = 0;
        yuv(x,y+1,0) =  yuv(x,y+1,1) = yuv(x,y+1,2) = 0;
        yuv(x+1,y+1,0) =  yuv(x+1,y+1,1) = yuv(x+1,y+1,2) = 0;

        output.push_back(gray);
        output.push_back(yuv);
    // Return the vector
    return output;
}


// White balances an image using the gray world assumption
FloatImage grayworld(const FloatImage & im)
{
    // Your code goes here
    FloatImage output(im);  // CHANGEME
    float sumr = 0, sumg = 0, sumb = 0;
    float ratio1;   // sumg/sumr
    float ratio2;   // sumg/sumb
    for (int x = 0; x < output.width(); ++x)
        for (int y = 0; y < output.height(); ++y)
        {   
            sumr = sumr + output(x,y,0);
            sumg = sumg + output(x,y,1);
            sumb = sumb + output(x,y,2);
        }  
    ratio1 = sumg/sumr;
    ratio2 = sumg/sumb;
    for (int x = 0; x < output.width(); ++x)
        for (int y = 0; y < output.height(); ++y)
        {  
            output(x,y,0) = output(x,y,0) * ratio1;
            output(x,y,1) = output(x,y,1) * ratio2;
        }
    return output;
}


// Histograms

// Stretches remaps the pixel values of channel c of im so that the minimum value maps to 0,
// and the maximum value maps to 1
void autoLevels(FloatImage & im, int c)
{
    // stretch pixel values of channel c to fill 0..1 range
    // you may find FloatImage::min and FloatImage::max useful
    float max = im.max();
    float min = im.min();
    float denominator = max - min;
    for (int x = 0; x < im.width(); ++x)
        for (int y = 0; y < im.height(); ++y)
        {  
            im(x,y,c) = ( im(x,y,c) - min ) / denominator;
        }
}


// initialize a histogram with numBins bins from the c-th channel of im
Histogram::Histogram(const FloatImage & im, int c, int numBins) :
    m_pdf(numBins, 0.0f), m_cdf(numBins, 0.0f)
{   
    // populate m_pdf with histogram values
    int index = 0;
    // float sum=0;
    m_cdf[0] = 0;   // Very important, otherwise, randomly allocated 
    FloatImage output(im);

    for (int x = 0; x < output.width(); ++x)
        for (int y = 0; y < output.height(); ++y)
        {     
             index = int(output(x,y,c) * (numBins - 1));
             m_pdf[index] = m_pdf[index] + 1;     // Careful about add interger to float
             if (index > numBins)
            cout<< index << '\n' << "Index Error in accessing m_pdf" <<endl;
        }

    for (int i = 0; i < numBins; i++)
     {  
         m_pdf[i] =  double(m_pdf[i] / (output.width() * output.height())) ;
        // sum += m_pdf[i];
        if (i == 0) 
            m_cdf[i] = m_pdf[i];
        else
            m_cdf[i] = m_cdf[i-1] + m_pdf[i];

        if (max_pdf < m_pdf[i])
            max_pdf = m_pdf[i]; 
 
     } 
     // cout<<"sum" <<sum<<endl;
     // cout <<"m_cdf = " << m_cdf[numBins-1] <<endl;
     // cout << "max_pdf =" << max_pdf << endl;
 
    // m_cdf = 
    // Grad/extra credit: populate m_cdf with running sum of m_pdf
}

// return the histogram bin index that value falls within
int Histogram::inverseCDF(double value) const
{
    int left = 0;
    int right = this->numBins();
    int index = (left + right) /2; 

    while (index != left && index != right) {   // Binary search 
        if ( value  > this->cdf(index) )
            left = index;                       // Left value preferable 
        else 
            right = index;
        index = (left + right) /2;
    }
    return index ; // CHANGEME
}


// Produce a numBins() x 100 x 3 image containing a visualization of the
// red, green and blue histogram pdfs
FloatImage visualizeRGBHistogram(const Histogram & histR,
                                 const Histogram & histG,
                                 const Histogram & histB)   // Find the max probility in three channels, caculate the desired pixels
{
    assert(histR.numBins() == histG.numBins() && histR.numBins() == histB.numBins());

    // create an image of appropriate size
    FloatImage im(histR.numBins(), 100, 3);
   // double base = max(histR.max_pdf,histG.max_pdf);
   //  base  = max(base,histB.max_pdf);
   //  double base = 0.1; 
   //  cout<< "base is" << base<<endl; 
   // cout<< histR.cdf(histR.numBins()-1) << " " <<histG.cdf(histR.numBins()-1)<<" "<<histB.cdf(histR.numBins()-1)<<endl;

    for (int i = 0; i < im.width(); ++i)    // use different bases for each channel; 
    {
        for(int j = 0; j < min(100 * histR.pdf(i) /histR.max_pdf, 100.0); ++j)  // Careful, the image is top down. The center point is on the top left.!!
            im(i,99 - j,0) = 1;
        
        for(int k = 0; k < min(100 * histG.pdf(i) /histG.max_pdf,100.0); ++k)
            im(i,99 - k,1) = 1;
        
        for(int l = 0; l < min(100 * histB.pdf(i) /histB.max_pdf,100.0); ++l)
            im(i,99 - l,2) = 1;

    }

    return im;
}


// Return a FloatImage which is the result of applying histogram equalization to im
FloatImage equalizeRGBHistograms(const FloatImage & im)
{
    int numLevels = 256;
    int index_r = 0 ;
    int index_g = 0 ;
    int index_b = 0 ;

    FloatImage output = im;
    Histogram His_r(im,0);
    Histogram His_g(im,1);
    Histogram His_b(im,2);

        for (int x = 0; x < im.width(); ++x)
            for (int y = 0; y < im.height(); ++y)
                 
        {   
            index_r = int(output(x,y,0) * (numLevels - 1));
            index_g = int(output(x,y,1) * (numLevels - 1));
            index_b = int(output(x,y,2) * (numLevels - 1)); 
            output(x,y,0) = (His_r.cdf(index_r) - His_r.pdf(0))/(1 - His_r.pdf(0));
            output(x,y,1) = (His_g.cdf(index_g) - His_g.pdf(0))/(1 - His_g.pdf(0));
            output(x,y,2) = (His_b.cdf(index_b) - His_b.pdf(0))/(1 - His_b.pdf(0));
       }
    // apply histogram equalization to each channel
    return output;
}

// Return a FloatImage which is the result of transfering the histogram of im2 onto the image data in im1
FloatImage matchRGBHistograms(const FloatImage & im1, const FloatImage & im2)   // Must imput the images with the same size.
{
    // perform histogram matching
    float numBins = 256.0;

    FloatImage after = im1;
    FloatImage before = im2;
    Histogram after_his_r(after,0);
    Histogram after_his_g(after,1);
    Histogram after_his_b(after,2);
    Histogram before_his_r(before,0);
    Histogram before_his_g(before,1);
    Histogram before_his_b(before,2);

    for (int x = 0; x < after.width(); ++x) // Careful!! Types conversion 
        for (int y = 0; y < after.height(); ++y){
            after(x,y,0) = float ( before_his_r.inverseCDF(after_his_r.cdf(int (after(x,y,0) * (numBins-1))))) / ((numBins - 1));  
            after(x,y,1) = float ( before_his_g.inverseCDF(after_his_g.cdf(int (after(x,y,1) * (numBins-1)))) / (numBins - 1) );
            after(x,y,2) = float ( before_his_b.inverseCDF(after_his_b.cdf(int (after(x,y,2) * (numBins-1)))) / (numBins - 1) );
            }

    // Histogram after_test(after,0);   // Test the processed image, calculate the differences of cdf. 
    //      for (int i = 0; i < 255; ++i)
    //  {
    //     cout << "difference" << after_test.cdf(i) - before_his_r.cdf(i) <<endl; 
    //  }

    return after;
}

