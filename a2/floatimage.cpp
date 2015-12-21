/*!
    \file floatimage.cpp
    \brief Contains the implementation of a floating-point image class with an arbitrary number of channels
    \author Wojciech Jarosz
    
    CS 89/189 Computational Aspects of Digital Photography C++ basecode.
*/

#include "floatimage.h"
#include "utils.h"
#include "lodepng.h"
#include <math.h>
#include <iostream>
#include <sstream>

using namespace std;

// anonymous (unnamed) namespaces are used in C++ to make functions and variable
// local to the file and not pollute the global namespace. This is the same as
// the C way of having a static global variable or static function
namespace
{

float uint8_to_float(const unsigned char in)
{
    return ((float) in)/(255.0f);
}

unsigned char float_to_uint8(float in)
{
    return (unsigned char) (255.0f*clamp(in, 0.0f, 1.0f));
}

void compareXYDimensions(const FloatImage & im1, const FloatImage & im2)
{
    if (im1.sizeX() != im2.sizeX() ||
        im1.sizeY() != im2.sizeY())
        throw MismatchedDimensionsException();
}

} // end namespace


FloatImage::FloatImage() : Array3D<float>()
{
    // empty
}


FloatImage::FloatImage(int width, int height, int channels) :
    Array3D<float>(width, height, channels)
{
    // empty
}


FloatImage::FloatImage(const string & filename) : Array3D<float>()
{
    read(filename);
}


FloatImage::FloatImage(const FloatImage &in)
{
    resize(in.width(), in.height(), in.depth());
    m_data = in.m_data;
}


FloatImage& FloatImage::operator=(const FloatImage &in)
{
    m_data = in.m_data;
    m_sizeX = in.m_sizeX;
    m_sizeY = in.m_sizeY;
    m_sizeZ = in.m_sizeZ;
    m_strideZ = in.m_strideZ;
    return *this;
}

void FloatImage::clear(const vector<float> & channelValues)
{
    // check z bounds

    for (int z = 0; z < sizeZ(); ++z)
        for (int y = 0; y < sizeY(); ++y)
            for (int x = 0; x < sizeX(); ++x)
                operator()(x,y,z) = channelValues[z];
}


FloatImage & FloatImage::operator+=(float s)
{
    for (int i = 0; i < size(); ++i)
        operator()(i) += s;
    return *this;
}


FloatImage & FloatImage::operator*=(float s)
{
    for (int i = 0; i < size(); ++i)
        operator()(i) *= s;
    return *this;
}


FloatImage & FloatImage::operator/=(float s)
{
    if (s == 0.0f)
        throw DivideByZeroException();
    return *this *= 1.0f/s;
}


FloatImage & FloatImage::operator+=(const FloatImage &other)
{
    compareXYDimensions(*this, other);

    if (this->channels() == other.channels())
    {
        for (int i = 0; i < size(); ++i)
            operator()(i) += other(i);
    }
    else
    {
        if (other.channels() != 1)
            throw MismatchedDimensionsException();

        // if other is a grayscale image, add it to each channel
        for (int c = 0; c < channels(); c++)
            for (int y = 0; y < height(); y++)
                for (int x = 0; x < width(); x++)
                    operator()(x,y,c) += other(x,y,0);
    }
    return *this;
}


FloatImage & FloatImage::operator-=(const FloatImage &other)
{
    compareXYDimensions(*this, other);

    if (this->channels() == other.channels())
    {
        for (int i = 0; i < size(); ++i)
            operator()(i) -= other(i);
    }
    else
    {
        if (other.channels() != 1)
            throw MismatchedDimensionsException();

        // if other is a grayscale image, subtract it from each channel
        for (int c = 0; c < channels(); c++)
            for (int y = 0; y < height(); y++)
                for (int x = 0; x < width(); x++)
                    operator()(x,y,c) -= other(x,y,0);
    }
    return *this;
}


FloatImage & FloatImage::operator*=(const FloatImage &other)    // Consider the dimension. 
{
    compareXYDimensions(*this, other);

    if (this->channels() == other.channels())
    {
        for (int i = 0; i < size(); ++i)
            operator()(i) *= other(i);
    }
    else
    {
        if (other.channels() != 1)          // Consider 1 channel for gray-scale picture 
            throw MismatchedDimensionsException();

        // if other is a grayscale image, multiply each channel by it
        for (int c = 0; c < channels(); c++)
            for (int y = 0; y < height(); y++)
                for (int x = 0; x < width(); x++)
                    operator()(x,y,c) *= other(x,y,0);
    }
    return *this;
}


FloatImage & FloatImage::operator/=(const FloatImage &other)
{
    compareXYDimensions(*this, other);

    if (this->channels() == other.channels())
    {
        for (int i = 0; i < size(); ++i)
        {
            if (other(i) == 0.0f)
                throw DivideByZeroException();
            operator()(i) /= other(i);
        }
    }
    else
    {
        if (other.channels() != 1)
            throw MismatchedDimensionsException();

        // if other is a grayscale image, divide each channel by it
        for (int c = 0; c < channels(); c++)
            for (int y = 0; y < height(); y++)
                for (int x = 0; x < width(); x++)
                {
                    if (other(x,y,0) == 0.0f)
                        throw DivideByZeroException();
                    operator()(x,y,c) /= other(x,y,0);
                }
    }
    return *this;
}

FloatImage operator*(float s, const FloatImage &other)
{
    FloatImage ret(other);
    return ret *= s;
}

FloatImage operator/(float s, const FloatImage &other)
{
    FloatImage ret(other.width(), other.height(), other.channels());
    for (int i = 0; i < ret.size(); ++i)
    {
        if (other(i) == 0.0f)
            throw DivideByZeroException();
        ret(i) = s/other(i);
    }
    return ret;
}

FloatImage operator+(float s, const FloatImage &other)
{
    FloatImage ret(other);
    return ret += s;
}

FloatImage operator-(float s, const FloatImage &other)
{
    FloatImage ret(other);
    for (int i = 0; i < ret.size(); ++i)
        ret(i) = s - ret(i);
    return ret;
}

FloatImage FloatImage::operator+(float s) const
{
    FloatImage ret(*this);
    return ret += s;
}

FloatImage FloatImage::operator*(float s) const
{
    FloatImage ret(*this);
    return ret *= s;
}

FloatImage FloatImage::operator+(const FloatImage &other) const
{
    FloatImage ret(*this);
    return ret += other;
}


FloatImage FloatImage::operator-(const FloatImage &other) const
{
    FloatImage ret(*this);
    return ret -= other;
}


FloatImage FloatImage::operator*(const FloatImage &other) const
{
    FloatImage ret(*this);
    return ret *= other;
}


FloatImage FloatImage::operator/(const FloatImage &other) const
{
    FloatImage ret(*this);
    return ret /= other;
}

float FloatImage::min(int c) const
{
    float mn = operator()(0,0,c);
    
    for (int y = 0; y < height(); y++)
        for (int x = 0; x < width(); x++)
            mn = std::min(mn, operator()(x,y,c));

    return mn;
}

float FloatImage::min() const
{
    float mn = 0;
    for (int c = 0; c < channels(); ++c)
        mn = std::min(mn, min(c));
    return mn;
}

float FloatImage::max(int c) const
{
    float mx = operator()(0,0,c);
    
    for (int y = 0; y < height(); y++)
        for (int x = 0; x < width(); x++)
            mx = std::max(mx, operator()(x,y,c));

    return mx;
}

float FloatImage::max() const
{
    float mx = 0;
    for (int c = 0; c < channels(); ++c)
        mx = std::max(mx, max(c));
    return mx;
}




bool FloatImage::read(const string & filename)
{
    vector<unsigned char> uint8_image;
    unsigned int height_;
    unsigned int width_;
    unsigned int channels_ = 4;
    unsigned int outputchannels_ = 3; // Throw away transparency

    // In column major order with packed color values
    unsigned error = lodepng::decode(uint8_image, width_, height_,
                                     filename.c_str());

    // if there's an error, display it
    if (error)
    {
        cerr << "PNG decoder error " << error << ": " << lodepng_error_text(error) << endl;
        return false;
    }

    resize(width_, height_, outputchannels_, 0.0f);

    for (unsigned int x = 0; x < width_; x++)
        for (unsigned int y = 0; y < height_; y++)
            for (unsigned int c = 0; c < outputchannels_; c++)
                operator()(x,y,c) = uint8_to_float(uint8_image[c + x*channels_ + y*channels_*width_]);  // TODO: Is this the same as that in mdata(a,b,c)?  array3d.h

    return true;
}


bool FloatImage::write(const string &filename) const
{
    if (channels() != 1 && channels() != 3 && channels() != 4)
        throw ChannelException();

    int png_channels = 4;
    vector<unsigned char> uint8_image(height()*width()*png_channels, 255);
    int c;
    for (int x = 0; x < width(); x++)
        for (int y = 0; y < height(); y++)
        {
            for (c = 0; c < channels(); c++)
                uint8_image[c + x*png_channels + y*png_channels*width()] = float_to_uint8(operator()(x,y,c));

            for ( ; c < 3; c++)
                // Only executes when there is one channel
                uint8_image[c + x*png_channels + y*png_channels*width()] = float_to_uint8(operator()(x,y,0));
        }

    unsigned error = lodepng::encode(filename.c_str(), uint8_image, width(), height());

    // if there's an error, display it
    if (error)
        cerr << "PNG encoder error " << error << ": "<< lodepng_error_text(error) << endl;
    return error == 0;
}

int FloatImage::s_debugWriteNumber = 0;

bool FloatImage::debugWrite() const
{
    ostringstream ss;
    ss << "./Output/" <<  s_debugWriteNumber << ".png";
    string filename = ss.str();
    s_debugWriteNumber++;
    return write(filename);
}


