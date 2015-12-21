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
#include "utils.h"
#include <math.h>

using namespace std;

// Basic denoising by computing the average of a sequence of images
FloatImage denoiseSeq(const vector<FloatImage> &imSeq)
{
    // create an appropriately sized output FloatImage
    FloatImage output(imSeq[0].width(), imSeq[0].height(), imSeq[0].channels());  //   CHANGEME
    int im_len = imSeq.size(); 
    float sum = 0.0; 

    if (im_len == 1)
        return imSeq[0];

    // Ensure they are perfectly alligned   TODO: HoW to ensure the alignment ? 
    // Ensure they have the same sizes
    for (int i = 1; i < im_len;i++)
        if (imSeq[i].width() != imSeq[i-1].width() || imSeq[i].height() != imSeq[i-1].height() || imSeq[i].channels() != imSeq[i-1].channels())
            throw MismatchedDimensionsException();

    //  compute average of sequence
    for (int n = 0; n < output.size(); n++){
        sum = 0;
        for (int j = 0; j < im_len; j++)
        sum += imSeq[j](n);
        output(n) =  sum / im_len; 
    }

     return output;
}


// returns an image visualizing the per-pixel and
// per-channel log of the signal-to-noise ratio scaled by scale.
FloatImage logSNR(const vector<FloatImage> &imSeq, float scale)
{
    // create an appropriately sized logSNR FloatImage
    FloatImage logSNR(imSeq[0].width(), imSeq[0].height(), imSeq[0].channels());
    int im_len = imSeq.size(); 
    float ave[logSNR.size()]; 
    float var[logSNR.size()]; 
    float lg = 0;


    if (im_len == 1)
        return imSeq[0];

    // Ensure they have the same sizes
    for (int i = 1; i < im_len;i++)
        if (imSeq[i].width() != imSeq[i-1].width() || imSeq[i].height() != imSeq[i-1].height() || imSeq[i].channels() != imSeq[i-1].channels())
            throw MismatchedDimensionsException();

    //  compute average of sequence
    for (int i = 0; i < logSNR.size(); i++){
        for (int j = 0; j < im_len; j++)
            ave[i] += imSeq[j](i);
         ave[i] = ave[i] / im_len ;
    }

    //  compute variance of sequence
    for (int i = 0; i < logSNR.size(); ++i){
        for (int j = 0; j < im_len; ++j){
             var[i] = var[i] + pow((ave[i] - imSeq[j](i)),2);}
            // Unbiased estimator
            var[i] = var[i] / (im_len - 1); 
          //   cout<< ave[i] <<"variance"<<var[i]<<endl;
            if (var[i] == 0.0)
              var[i] = 1e-6;
        // set max variance as numbers of images
            if (var[i] > im_len)
             var[i] =  30;
        //  computer 10 x log of the signal-to-noise ratio multiplied by scale   
        lg = 5 * scale * log10((ave[i] / var[i]));
                  //   cout<< ave[i] <<"variance"<<var[i]<<endl;
        logSNR(i) = lg;
        //cout<<logSNR(i)<<endl; 
            
    }
    
    return logSNR;
}


// returns the (x,y) offset that best aligns im2 to match im1.
vector<int> align(const FloatImage &im1, const FloatImage &im2, int maxOffset)
{
    //  offset that best aligns images
    vector<int> offset(2, 0);       // Don't let it be  0  0. !!!
    FloatImage img(im1.width(),im1.height(),im1.channels());
    float sum_sqr_dif = 0.0;
    float min_sum_sqr_dif = 0.0;

    //  find the best offset to align im1 to im2 for all the pixels <= maxOffset
    //  Hint: you might find roll() useful
    for (int x = - maxOffset; x <= maxOffset; ++x)
        for (int y = - maxOffset; y <= maxOffset; ++y){

            sum_sqr_dif = 0.0;
            img = im1;
            img = roll(img,x,y);  // Don't roll im1!!!!
            for (int i = maxOffset + 1; i < img.width() - maxOffset; ++i)
                for (int j = maxOffset + 1; j < img.height() - maxOffset; ++j)
                    for(int k = 0; k < img.channels();++k){
                    sum_sqr_dif = sum_sqr_dif + pow((im2(i,j,k) - img(i,j,k)),2);
                    // cout<<"x="<<x<<"y="<<y<<"sum_sqr_dif="<<sum_sqr_dif<<"min_sum_sqr_dif="<<min_sum_sqr_dif<<endl;
                }
            if (min_sum_sqr_dif > sum_sqr_dif){
                    min_sum_sqr_dif = sum_sqr_dif;
                offset[0] = x;
                offset[1] = y;
                }
             if ( (x == - maxOffset) && (y == - maxOffset)){     
            // Wired!!  Put it in of the Line 130, It returns "illegal instruction:4"; Compiler Optimization bugs???
            // TODO: apply command like  mmacosx-version-min=10.5
                min_sum_sqr_dif = sum_sqr_dif;
            } 

            }
           // cout << "sum_sqr is" << sum_sqr_dif <<"min is"<<min_sum_sqr_dif <<"(x,y) is"<<x<<','<<y<<endl;    
    return offset;
}

// registers all images to the first one in a sequence and outputs
// a denoised image even when the input sequence is not perfectly aligned.
FloatImage alignAndDenoise(const vector<FloatImage> &imSeq, int maxOffset)
{
    //  register all images to the first one in sequence
    FloatImage output(imSeq[0].width(),imSeq[0].height(),imSeq[0].channels());  
    vector<FloatImage> imSeq_copy(imSeq);
    int im_len = imSeq.size(); 
    vector<int> offset(2,0);

    //  iterate over the remaining images in sequence
    //  and do the required processing to denoise image by aligning with first image in sequence
    for (int i = 1; i < im_len; ++i){
        offset = align(imSeq_copy[i],imSeq_copy[0],maxOffset);
        imSeq_copy[i] = roll(imSeq_copy[i],offset[0],offset[1]);
        // cout<<offset[0]<<'-'<<offset[1]<<endl;
    }

    output = denoiseSeq(imSeq_copy);
    return output;
}

// split a Sergey images to turn it into one 3-channel image.
FloatImage split(const FloatImage &sergeyImg)
{
    //  create output image with appropriate height (hint: use floor())
    FloatImage output(sergeyImg.width(),sergeyImg.height() / 3,3);  //  CHANGEME
    FloatImage sergeyImg_copy(sergeyImg);
    int height_1 = output.height();
    int height_2 = 2 * height_1;

    //  store values from sergeyImg into output's 3 channels
    for (int i=0; i< output.width();++i){
        for (int j = 0; j < height_1 ; ++j){
            output(i,j,0) = sergeyImg_copy(i,j,0);
            output(i,j,1) = sergeyImg_copy(i,j + height_1,1);
            output(i,j,2) = sergeyImg_copy(i,j + height_2,2);
    }
    }


    return output;
}

// aligns the green and blue channels of your rgb channel of a sergey
// image to the red channel. This should return an aligned RGB image
FloatImage sergeyRGB(const FloatImage &sergeyImg, int maxOffset)
{
    //  call split()
    FloatImage rgb(sergeyImg); //  CHANGEME
    FloatImage img(rgb);
    FloatImage output(sergeyImg.width(),sergeyImg.height()/3,3);
    FloatImage green_channel(sergeyImg);
    FloatImage blue_channel(sergeyImg);

    vector<int> g_position(2,0);
    vector<int> b_position(2,0);
    float sum_sqr_dif_g;
    float min_sum_sqr_dif_g;
    float min_sum_sqr_dif_b;
    float sum_sqr_dif_b;

    //  align green to red channel 
    //  align blue to red channel

    for (int x = - maxOffset; x <= maxOffset; ++x)
        for (int y = - maxOffset; y <= maxOffset; ++y){
            img = rgb; 
            sum_sqr_dif_g = 0.0;
            sum_sqr_dif_b = 0.0;
            img = roll(img,x,y);
            

            for (int i = maxOffset + 1; i < rgb.width() - maxOffset; ++i)
                for (int j = maxOffset +1; j < (rgb.height() - maxOffset )/ 3; ++j){
                    sum_sqr_dif_g = sum_sqr_dif_g + pow((rgb(i,j,0) - img(i,output.height() + j,1)),2);
                    sum_sqr_dif_b = sum_sqr_dif_b + pow((rgb(i,j,0) - img(i,2 * output.height() + j,2)),2);
                }
            cout<<"min_sum_sqr_dif_g"<<min_sum_sqr_dif_g<<" min_sum_sqr_dif_b"<<min_sum_sqr_dif_b<<endl;
                    
            if (min_sum_sqr_dif_g > sum_sqr_dif_g){
                min_sum_sqr_dif_g = sum_sqr_dif_g;
                g_position[0] = x;
                g_position[1] = y;
            }
            if (min_sum_sqr_dif_b > sum_sqr_dif_b){
                min_sum_sqr_dif_b = sum_sqr_dif_b;
                b_position[0] = x;
                b_position[1] = y;
            }
            if ( (x == - maxOffset) && (y == - maxOffset+1)){     
            // Wired!!  Put it in of the Line 130, It returns "illegal instruction:4"; Compiler Optimization bugs???
            // TODO: apply command like  mmacosx-version-min=10.5
                min_sum_sqr_dif_g = sum_sqr_dif_g;
                min_sum_sqr_dif_b = sum_sqr_dif_b;
                cout<< "yep"<<min_sum_sqr_dif_g<<endl;
                cout<< "yep"<<min_sum_sqr_dif_b<<endl;

            }

         }

         green_channel = roll(green_channel,g_position[0],g_position[1]);
         blue_channel = roll(blue_channel,b_position[0],b_position[1]);
         cout<<"ok here"<<endl;
         cout<<"gx"<<g_position[0]<<"gy"<<g_position[1]<<endl;
         cout<<"bx"<<b_position[0]<<"by"<<b_position[1]<<endl;
    for (int x = 0; x < output.width() ; ++x)
        for (int y = 0; y < output.height() ; ++y){
            output(x,y,0) = rgb(x,y,0);
            output(x,y,1) = green_channel(x,output.height() + y,1);
            output(x,y,2) = blue_channel(x,2 * output.height() + y,2);
        }


    return output;
}


/**************************************************************
 //               DON'T EDIT BELOW THIS LINE                //
 *************************************************************/

// This circularly shifts an image by xRoll in the x direction and
// yRoll in the y direction. xRoll and yRoll can be negative or postive
FloatImage roll(const FloatImage &im, int xRoll, int yRoll)
{    
    int xNew, yNew;
    FloatImage imRoll(im.width(), im.height(), im.channels());
    
    // for each pixel in the original image find where its corresponding
    // location is in the rolled image
    for (int x = 0; x < im.width(); x++)
    {
        for (int y = 0; y < im.height(); y++)
        {
            // use modulo to figure out where the new location is in the
            // rolled image. Then take care of when this returns a negative number
            xNew = (x + xRoll) % im.width();
            yNew = (y + yRoll) % im.height();
            xNew = (xNew<0)*(imRoll.width() + xNew) + (xNew>=0)*xNew;
            yNew = (yNew<0)*(imRoll.height() + yNew) + (yNew>=0)*yNew;
            
            // assign the rgb values for each pixel in the original image to
            // the location in the new image
            for (int z = 0; z < im.channels(); z++)
                imRoll(xNew, yNew, z) = im(x,y,z);
        }
    }
    
    // return the rolled image
    return imRoll;
}
