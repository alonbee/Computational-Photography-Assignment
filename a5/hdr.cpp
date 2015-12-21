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

// hdr.cpp
// Assignment 5


#include "hdr.h"
#include "filtering.h"
#include <math.h>
#include <algorithm>


using namespace std;

/**************************************************************
 //                       HDR MERGING                        //
 *************************************************************/

// Generate a weight image that indicates which pixels are good to use in hdr
FloatImage computeWeight(const FloatImage &im, float epsilonMini, float epsilonMaxi){
  FloatImage mask(im.width(),im.height(),im.channels());
  for (int i =0; i < mask.width(); ++i)
    for (int j = 0; j < mask.height(); ++j)
      for (int k = 0; k < mask.channels(); ++k)
      {
      mask(i,j,k) = im(i,j,k) > epsilonMaxi ? 0.0 : im(i,j,k) < epsilonMini ? 0.0 : 1.0; 
      }

    return mask; 
}

// Compute the multiplication factor between a pair of images
float computeFactor(const FloatImage &im1, const FloatImage &w1, const FloatImage &im2, const FloatImage &w2){
  FloatImage m1(im1);
  FloatImage msk1(w1);
  FloatImage m2(im2);
  FloatImage msk2(w2);
  float median = 0;
  float value = 0;
  std::vector<float> val;

  for (int i = 0; i < m1.width(); ++i)
    for (int j = 0; j < m1.height(); ++j)
      for (int k = 0; k < m1.channels(); ++k)
         {
           if (msk1(i,j,k) == 1 && msk2(i,j,k) == 1)
           {
            if (m1(i,j,k) < 10e-10)    
              m1(i,j,k) = 10e-10;
            value = (float) (m2(i,j,k) / m1 (i,j,k));
            val.push_back(value);
           }
         }
    std::sort (val.begin(), val.end());             
    median = val[val.size() / 2]; 
    return median;  
}

// Merge images to make a single hdr image
FloatImage makeHDR(vector<FloatImage> &imSeq, float epsilonMini, float epsilonMaxi){

    // invert gamma correction
    for (int i = 0; i < ((int) imSeq.size()) ; i++) {
      imSeq[i] = changeGamma(imSeq[i], 1.0/2.2, 1.0f);
    }
  std::vector<FloatImage> images(imSeq);
  std::vector<FloatImage> masks(imSeq);
  FloatImage output(images[0]);
  int max = imSeq.size() - 1;
  int min = 0; //TODO: BUG?

  float sum = 0;
  float sum_mask = 0;
  float coef[imSeq.size()];

  for (int m = 0; m < ((int) imSeq.size()) ; ++m){
      masks[m] =  computeWeight(images[m],epsilonMini,epsilonMaxi);
  }
// special cases for darkest and brightest images
  for (int i = 0; i < images[min].size();++i){
    if (images[min](i) >= epsilonMaxi)    // TODO: threshold in one direction for these cases
      masks[min](i) = 1.0;
    if (images[max](i) <= epsilonMini)
      masks[max](i) = 1.0;
  }

  coef[0] = 1;
  for (int i = 1; i < imSeq.size() ; ++i){
    // images[i-1] = images[i-1] + 10e-10;  // Quite important to avoid zero in calculating luminance and chrominance
    coef[i] = computeFactor(images[i-1] ,masks[i-1],images[i],masks[i]);
  }
     cout<<"HDR step 1 finished" <<endl;

  // images[imSeq.size()-1] = images[imSeq.size()-1] + 10e-10;
  // cout<<"len"<<images.size()<<endl;
  // cout<<"coef[0] "<<coef[0]<<endl;

  for (int i = imSeq.size()-1;i >0; --i){   
        for (int j = i-1; j >= 0;  --j)
          coef[i] = coef[i] * coef[j];
        // cout<<"coef["<<i<<"] "<<coef[i]<<endl;
  }
       cout<<"HDR step 2 finished" <<endl;

  for(int i = 0; i < images[0].width(); ++i)
    for(int j = 0; j < images[0].height(); ++j)
      for(int k = 0; k < images[0].channels(); ++k){
        sum = 0;
        sum_mask = 0;
          for(int m = 1; m < ((int) imSeq.size()) ; ++m)
        {
          // if (coef [m] == 0)
          //   coef [m] = 10e-10;
          sum += float(masks[m](i,j,k) * images[m](i,j,k) / coef[m]);
          sum_mask += masks[m](i,j,k);
        }
        if (sum_mask == 0)
          output(i,j,k) = images[0](i,j,k);
        else
          output(i,j,k) = ((images[0](i,j,k) * masks[0](i,j,k) / coef[0]) + sum) / (sum_mask + masks[0](i,j,k));
      }
    //output = changeGamma(output, 1.0f, 1.0/2.2);

    return output;  
}

/**************************************************************
 //                      TONE MAPPING                        //
 *************************************************************/


// Tone map an hdr image
FloatImage toneMap(const FloatImage &im, float targetBase, float detailAmp, bool useBila, float sigmaRange) {

  // UNCOMMENT THIS LINE!
  // add gamma correction back into the image right before returning
  //outImage = changeGamma(outImage, 1.0f, 1/2.2);
  float min_lumi; 
  float sigmaDomain;
  float targ;
  float range;
  float factor;
  float ma;
  float mi;
  FloatImage gaussian_bilater(im);


  FloatImage output(im.width(),im.height(),im.channels());
  std::vector<FloatImage> lumi_chrom(2);

  lumi_chrom = lumiChromi(im);
 
  min_lumi= image_minnonzero(lumi_chrom[0]);

  FloatImage lumi_log(lumi_chrom[0] + min_lumi) ;
  lumi_log = log10FloatImage(lumi_log);
    
  FloatImage lumi_log_output(lumi_log.width(),lumi_log.height(),lumi_log.channels());
  FloatImage detail(lumi_log.width(),lumi_log.height(),lumi_log.channels());
  FloatImage im_after(lumi_log.width(),lumi_log.height(),lumi_log.channels());
  sigmaDomain = (float) max(im.height(),max(im.width(),im.channels())) / 50;
    // sigmaDomain = (float) max(im.height(),max(im.width(),im.channels())) / 500; 
  cout<<"step 1 finished" <<endl;

  if (useBila)
    gaussian_bilater = bilateral(lumi_log,sigmaRange,sigmaDomain,3,true);
  else
    gaussian_bilater = gaussianBlur_2D(lumi_log,sigmaDomain,3,true);

  cout<<"step 2 finished" <<endl;

  detail = lumi_log - gaussian_bilater;
  // cout<<"gaussian_bilater max = "<<gaussian_bilater.max()<<"min "<<gaussian_bilater.min()<<endl;

  targ = log10(targetBase);
  ma = gaussian_bilater.max();
  mi = gaussian_bilater.min();
  range = ma - mi;
  factor = (float) (targ / range); 

  // if (lumi_log.max() > 0)    // max base value mapped to 1 
  //   lumi_log = lumi_log - lumi_log.max();

  detail = detail * detailAmp;
  gaussian_bilater = factor * (gaussian_bilater - gaussian_bilater.max());
  im_after = gaussian_bilater + detail;
    cout<<"step 3 finished" <<endl;

  // im_after = lumi_log + detail;
  for (int i=0; i< im_after.size();i++)
    if (im_after(i) > 0)
      im_after(i) = 0;

  // Caculate the offset.
  
  lumi_log_output = exp10FloatImage(im_after);
  output = lumiChromi2rgb(lumi_log_output,lumi_chrom[1]);

  output = changeGamma(output, 1.0f, 1/2.2);

  return output; // change this
}



// Tone Mapping Helpful Functions

// image --> log10FloatImage
FloatImage log10FloatImage(const FloatImage &im) {
  // Taking a linear image im, transform to log10 scale.
  // To avoid infinity issues, make any 0-valued pixel be equal the the minimum
  // non-zero value. See image_minnonzero(im).
  FloatImage img(im);
  FloatImage output(im.width(),im.height(),im.channels());

  float min_val = image_minnonzero(img);

  for (int i = 0; i < im.size(); ++i){
    if (img(i) == 0)
      img(i) = min_val;
    output(i) =  log10(img(i)); 
  }
  return output; // change this
}

// FloatImage --> 10^FloatImage
FloatImage exp10FloatImage(const FloatImage &im) {
  // take an image in log10 domain and transform it back to linear domain.
  // see pow(a, b)
  FloatImage output(im.width(),im.height(),im.channels());
  for (int i = 0; i < im.size(); ++i)
    output(i) = pow(10,im(i));
  return output; // change this
}

// min non-zero pixel value of image
float image_minnonzero(const FloatImage &im) {
  float min = 10e+30;
  for (int i = 0; i < im.size(); ++i)
    if(im(i) < min && im(i) != 0)
      min = im(i);
    // cout<< "min = "<<min<<endl;
  return min; // change this
}
