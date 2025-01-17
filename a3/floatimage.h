/*!
    \file floatimage.h
    \brief Contains the definition of a floating-point image class with an arbitrary number of channels
    \author Wojciech Jarosz
    
    CS 89/189 Computational Aspects of Digital Photography C++ basecode.
*/
#ifndef __FLOATIMAGE_H
#define __FLOATIMAGE_H

#include "array3D.h"
#include <string>


//! Floating point image
class FloatImage : public Array3D<float>
{
public:
    //-----------------------------------------------------------------------
    //@{ \name Constructors, destructors, etc.
    //-----------------------------------------------------------------------
    FloatImage();
    FloatImage(int width, int height, int channels);
    FloatImage(const std::string & filename);
    FloatImage(const FloatImage &other);
    FloatImage & operator=(const FloatImage &other);
    //@}


    //-----------------------------------------------------------------------
    //@{ \name Pixel accessors & manipulators.
    //-----------------------------------------------------------------------
    int channels() const {return m_sizeZ;}
    void clear(float v = 0.0f) {reset(v);}
    void clear(const std::vector<float> & channelValues);
    //@}


    //-----------------------------------------------------------------------
    //@{ \name Component-wise arithmetic and assignment.
    //-----------------------------------------------------------------------
    friend FloatImage operator*(float s, const FloatImage &other);
    friend FloatImage operator/(float s, const FloatImage &other);
    friend FloatImage operator+(float s, const FloatImage &other);
    friend FloatImage operator-(float s, const FloatImage &other);
    
    FloatImage operator*(float s) const;
    FloatImage operator/(float s) const {return this->operator*(1.0f/s);}
    FloatImage operator+(float s) const;
    FloatImage operator-(float s) const {return this->operator+(-s);}

    FloatImage operator+(const FloatImage &other) const;
    FloatImage operator-(const FloatImage &other) const;
    FloatImage operator*(const FloatImage &other) const;
    FloatImage operator/(const FloatImage &other) const;

    FloatImage & operator+=(float s);
    FloatImage & operator-=(float s) {return *this += -s;}
    FloatImage & operator*=(float s);
    FloatImage & operator/=(float s);

    FloatImage & operator+=(const FloatImage &other);
    FloatImage & operator-=(const FloatImage &other);
    FloatImage & operator*=(const FloatImage &other);
    FloatImage & operator/=(const FloatImage &other);
    //@}


    //-----------------------------------------------------------------------
    //@{ \name Image statistics.
    //-----------------------------------------------------------------------
    float min() const;
    float min(int channel) const;
    float max() const;
    float max(int channel) const;
    //@}


    //-----------------------------------------------------------------------
    //@{ \name Image I/O.
    //-----------------------------------------------------------------------
    bool read(const std::string &name);
    bool write(const std::string &name) const;
    bool debugWrite() const; // Writes image to Output directory with automatically chosen name
    //@}

private:
    static int s_debugWriteNumber; // Image number for debug write
};

#endif
