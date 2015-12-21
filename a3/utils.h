/*!
    \file utils.h
    \brief Contains the definition of various utility functions.
    \author Wojciech Jarosz
    
    CS 89/189 Computational Aspects of Digital Photography C++ basecode.
*/
#ifndef __UTILS_H
#define __UTILS_H

//! Linear interpolation.
/*!
    Linearly interpolates between \a a and \a b, using parameter t.
    \param a A value.
    \param b Another value.
    \param t A blending factor of \a a and \a b.
    \return Linear interpolation of \a a and \b -
            a value between a and b if \a t is between 0 and 1.
*/
template <typename T, typename S>
inline T
lerp(T a, T b, S t)
{
    return T((S(1)-t) * a + t * b);
}


//! Inverse linear interpolation.
/*!
    Computes the interpolation factor t such that lerp(a,b,t) = m
*/
template <typename T>
inline T
lerpFactor(T a, T b, T m)
{
    return (m - a) / (b - a);
}


//! Clamps a value between two bounds.
/*!
    \param a The value to test.
    \param l The lower bound.
    \param h The upper bound.
    \return The value \a a clamped to the lower and upper bounds.
    
    This function has been specially crafted to prevent NaNs from propagating.
*/
template <typename T>
inline T
clamp(T a, T l, T h)
{
    return (a >= l) ? ((a <= h) ? a : h) : l;
}


//! Returns a modulus b.
template <typename T>
inline T
mod(T a, T b)
{
    int n = (int)(a/b);
    a -= n*b;
    if (a < 0)
        a += b;
    return a;
}


//! Returns a modulus 1, assuming a is positive.
template <typename T>
inline T
mod1(T a)
{
    return a - int(a);
}


template <typename T> inline T pow2 (T x) {return x*x;}
template <typename T> inline T pow3 (T x) {return x*x*x;}
template <typename T> inline T pow4 (T x) {T x2 = x*x; return x2*x2;}
template <typename T> inline T pow5 (T x) {T x2 = x*x; return x2*x2*x;}
template <typename T> inline T sqr  (T x) {return pow2(x);}
template <typename T> inline T cube (T x) {return pow3(x);}



#endif  // __UTILS_H
