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

// morphing.cpp
// Assignment 6


#include "morphing.h"
using namespace std;

/**************************************************************
 //            IMAGE WARPING/MORPHING FUNCTIONS              //
 *************************************************************/

// warp an entire image according to a pair of segments.
FloatImage warpBy1(const FloatImage &im, Segment &segBefore, Segment &segAfter)
{
    FloatImage output(im.width(),im.height(),im.channels());
    std::vector<float> X_before(2,0);
    float u;
    float v;

    for (int i = 0; i < output.width(); ++i)
        for(int j = 0; j < output.height();++j)
            for(int k = 0; k < output.channels();++k){
                u = segAfter.getU(i,j);    // From before or after??
                v = segAfter.getV(i,j);
                X_before = segBefore.UVtoX(u,v);
                output(i,j,k) = interpolateLin(im,X_before[0],X_before[1],k,true);
           //     output(i,j,k) = im.smartAccessor(X_before[0],X_before[1],k,true);
                // if (output(i,j,k)>=1)
                //     cout << "i ,j k"<<i<<","<<j<<","<<k<<"x="<<X_before[0]<<"y="<<X_before[1]<<endl;
            }
    return output; 
}

// warp an image according to a vector of before and after segments using segment weighting
FloatImage warp(const FloatImage &im, vector<Segment> &segsBefore, vector<Segment> &segsAfter, float a, float b, float p)
{   
    FloatImage output(im.width(),im.height(),im.channels());
    int size = segsBefore.size();
    std::vector<float> X_before(2,0);
    std::vector<float> X_after(2,0);
    float u;
    float v;
    float val_m;
    float weight;
    float weightsum;

    for (int i = 0; i < output.width(); ++i)
        for ( int j = 0; j < output.height(); ++j)
            for (int k = 0; k < output.channels(); ++k){
                val_m = 0;
                weightsum = 0;
                X_after[0] = 0;
                X_after[1] = 0;
                for (int m = 0; m < size; ++m){
                    u = segsAfter[m].getU(i,j);
                    v = segsAfter[m].getV(i,j);
                    X_before = segsBefore[m].UVtoX(u,v);
                    weight = segsAfter[m].weight(i,j,a,b,p);
                    X_after[0] += X_before[0] * weight;
                    X_after[1] += X_before[1] * weight;
                   // val_m += weight * interpolateLin(im,X_before[0],X_before[1],k,true);
                    // val_m += weight * im.smartAccessor(X_before[0],X_before[1],k,true);
                    weightsum += weight;
                }
                output(i,j,k) = interpolateLin(im,X_after[0]/weightsum,X_after[1]/weightsum,k,true);
            }
                

    return output; // CHANGE ME
}

// return a vector of N images in addition to the two inputs that morphs going between im1 and im2 for the corresponding segments
vector<FloatImage> morph(const FloatImage &im1, const FloatImage &im2, vector<Segment> &segsBefore, vector<Segment> &segsAfter, int N, float a, float b, float p)
{   
    if (im1.size() != im2.size() || segsBefore.size() != segsAfter.size())
        throw MismatchedDimensionsException();
    // check if im1 and im2 are the same size and throw the appropriate exception before performing morph
    int seg_size = segsBefore.size();
    float t; 
    FloatImage im_ini(im1.width(),im1.height(),im1.channels());
    vector<FloatImage> series(N+2,im_ini);
    vector<float> inter_vec;
    Segment sg(0.0,0.0,0.0,0);
    vector<Segment> sges_t(seg_size,sg);
    float inter_x1;
    float inter_y1;
    float inter_x2;
    float inter_y2;

    series[0] = im1;
    series[N+1] = im2;

    for (int l = 1; l < N+1; ++l){
     t = float(l) / (N + 1);
        cout<< "t = " <<t<<endl;
        FloatImage back_now_img(im1.width(),im1.height(),im1.channels());
        FloatImage future_now_img(im1.width(),im1.height(),im1.channels());
        FloatImage inter_img(im1.width(),im1.height(),im1.channels());

        for(int k = 0; k < seg_size;++k){
            inter_x1 = segsAfter[k].P[0] * t + segsBefore[k].P[0] * (1-t);

            inter_y1 = segsAfter[k].P[1] * t + segsBefore[k].P[1] * (1-t);
            inter_x2 = segsAfter[k].Q[0] * t + segsBefore[k].Q[0] * (1-t);
            inter_y2 = segsAfter[k].Q[1] * t + segsBefore[k].Q[1] * (1-t);
            // cout<<inter_x1<<""<<inter_y1<<""<<inter_x1<<""<<inter_x2<<endl;
            Segment intermidiary(inter_x1,inter_y1,inter_x2,inter_y2);
            sges_t[k] = intermidiary;
            // cout<< "sges_t[k] pq = ("<< sges_t[k].PQ[0]<<","<<sges_t[k].PQ[1]<<")"<<endl;
        }

        back_now_img = warp(im1,segsBefore,sges_t,a,b,p);
        future_now_img = warp(im2,segsAfter,sges_t,a,b,p);
        cout<<"ok"<<endl;

        inter_img = t * future_now_img + (1-t) * back_now_img;
        series[l] = inter_img;  

    } 
    return series; // CHANGE ME
}

/**************************************************************
 //                 SEGMENT CLASS FUNCTIONS                  //
 *************************************************************/

// Segment constructor takes in 2 points (x1,y1) and (x2,y2) correspoding to the ends of a segment and computes:
// P - 2-element vector to point (x1, y1)
// Q - 2-element vector to pont (x2, y2)
// PQ - 2-element vector from P to Q
// PQ2 - squared distance between P and Q
// PQlength - distance between P and Q
// PQDivByPQlength2 - 2-element vector PQ normalized by PQ2
// perpPQ - 2-element vector perpendicular to PQ
// perpPQDivByPQlength - 2-element vector perpPQ normalized by PQlength
Segment::Segment(float x1, float y1, float x2, float y2)
{
    P       = vector<float>(2,0);
    Q       = vector<float>(2,0);
    PQ      = vector<float>(2,0);
    perpPQ  = vector<float>(2,0);
    PQDivByPQlength2 = vector<float>(2,0);
    perpPQDivByPQlength = vector<float>(2,0);

    P[0] = x1;
    P[1] = y1;
    Q[0] = x2;
    Q[1] = y2;
    PQ = subtract(P,Q); 
    perpPQ[0] = - PQ[1];
    perpPQ[1] = PQ[0];
    PQ2 = dot(PQ,PQ);
    PQlength = sqrt(PQ2);
    PQDivByPQlength2 = scalarMult(PQ, 1.0 / PQ2);
    perpPQDivByPQlength = scalarMult(perpPQ, 1.0 / PQlength);
}


// Implement the computation of the u coordinate of a point (x,y) with respect to a segment
float Segment::getU(float x, float y)
{
    std::vector<float> PX(2,0);
    PX[0] = x - this->P[0];
    PX[1] = y - this->P[1];
    float u = 0.0;

    u = dot(PX,this->PQDivByPQlength2);
    return u; // CHANGE ME
}


// Implement the computation of the v coordinate of a point (x,y) with respect to a segment
float Segment::getV(float x, float y)
{
    std::vector<float> PX(2,0);
    PX[0] = x - this->P[0];
    PX[1] = y - this->P[1];
    float v = 0.0;
    v = dot(PX,this->perpPQDivByPQlength);
    return v; // CHANGE ME
}


// compute the new (x, y) position of a point given by the (u,v) location relative to another segment.
// return the point (x,y) in a 2-element vector
vector<float> Segment::UVtoX(float u, float v)
{
//takes the u,v values and returns a coordinate - to be called from target segment
    std::vector<float> x_new(2,0);
    x_new[0] = P[0] + u * PQ[0] + v * this->perpPQDivByPQlength[0];
    x_new[1] = P[1] + u * PQ[1] + v * this->perpPQDivByPQlength[1];
    // if (x_new[0] < 0 || x_new[1] <0 )
    //     cout<<"Negative values exit"<<endl;
    return x_new;  
}

// Implement distance from a point (x,y) to the segment. Remember the 3 cases from class
float Segment::dist(float x, float y)
{
    
    std::vector<float> PX(2,0);
    std::vector<float> QX(2,0);
    float distance;
    PX[0] = x - this->P[0];
    PX[1] = y - this->P[1];

    QX[0] = x - this->Q[0];
    QX[1] = y - this->Q[1];

    if (dot(this->PQ,PX) >= 0 && dot(this->PQ,QX) <= 0)
        distance = fabs(this->getV(x,y));
    else if(dot(PQ,PX) < 0)
        distance = sqrt((x - P[0]) * (x - P[0]) + (y - P[1]) * (y - P[1])); 
    else if (dot(PQ,QX) > 0)
        distance = sqrt((x - Q[0]) * (x - Q[0]) + (y - Q[1]) * (y - Q[1])); 
    else
        distance = 0;

    return distance; 
}


// compute the weight of a segment to a point (x,y) given the weight parameters a,b, and p
float Segment::weight(float x, float y, float a, float b, float p)
{
    float weight;
    weight = pow((float) pow(this->PQlength,p) / (a + this->dist(x,y)),b);   
    return weight; // CHANGE ME
}

/**************************************************************
 //               DON'T EDIT BELOW THIS LINE                //
 *************************************************************/

// subtracts 2 vectors of the same length.
vector<float> Segment::subtract(const vector<float> &vec1, const vector<float> &vec2)
{
// create vector from vec1 to vec2

    vector<float> vec_1_to_2 (vec1.size(), 0);

    if(vec1.size() == vec2.size()){

        for (int i=0; i<vec1.size(); i++){
            vec_1_to_2[i] = vec2[i] - vec1[i];
        }

    }else{
        throw MismatchedSizeException();
    }

    return vec_1_to_2;
}

// takes the dot product between 2 vectors of the same length
float Segment::dot(const vector<float> &vec1, const vector<float> &vec2)
{

    float dotProd = 0;

    if(vec1.size() == vec2.size()){

        for (int i=0; i<vec1.size(); i++){
            dotProd += vec2[i]*vec1[i];
        }

    }else{
        throw MismatchedSizeException();
    }

    return dotProd;
}

// mutliplies an entire vector by a scalor value
vector<float> Segment::scalarMult(const vector<float> &vec, float factor)
{

    vector<float> nVec (vec.size(), 0);
    for(int i=0; i<vec.size(); i++){
        nVec[i] = vec[i]*factor;
    }
    return nVec;
}

// destructor
Segment::~Segment() { } // Nothing to clean up
