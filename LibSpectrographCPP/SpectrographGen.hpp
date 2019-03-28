//
//  SpectrographGen.hpp
//  LibSpectrographCPP
//
//  Created by Mayo Furgerson on 3/27/19.
//  Copyright © 2019 WimMa Games. All rights reserved.
//

#ifndef SpectrographGen_hpp
#define SpectrographGen_hpp
#include <SFML/Audio.hpp>
#include <SFML/Graphics.hpp>
#include <complex>
#include <vector>
#include <math.h>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>

class SpectrographGen {
public:
    SpectrographGen();
    SpectrographGen(int SpecSampleRate,int SpecResPerSample);
    
    bool loadFromFile(std::string File);
    
    void saveToFile(std::string File);
protected:
    sf::Image generateSpectogram();
    void sort(std::complex<float>* a,int n);
    void fft(std::complex<float>* x,int n,int s);
    
    std::vector<sf::Int16> m_Samples;
    sf::Int32 m_SampleRate;
    
    int m_SpecSampleRate;
    int m_SpecResPerSample;
    
    int m_FFTSize;
    std::complex<float> e,t;
    std::vector<float> m_Hanning;
    std::vector<std::complex<float>> m_Temp;
    std::vector<std::vector<std::complex<float>>> m_Polar;
};

class SpectrographDecode {
public:
    SpectrographDecode();
    SpectrographDecode(int SpecSampleRate,int SampleRate,int ChannelCount);
    
    bool loadFromFile(std::string File);
    
    void saveToFile(std::string File);
protected:
    sf::SoundBuffer generateBuffer();
    void sort(std::complex<float>* a,int n);
    void ifft(std::vector<std::complex<float>>& x,int s);
    void fft(std::complex<float>* x,int n,int s);
    
    sf::Image m_Image;
    sf::Int32 m_SampleRate;
    int m_SpecSampleRate;
    int m_SpecResPerSample;
    int m_ChannelCount;
    
    int m_FFTSize;
    std::complex<float> e,t;
    std::vector<float> m_Hanning;
    std::vector<std::complex<float>> m_Temp;
    std::vector<std::vector<std::complex<float>>> m_Polar;
};




/*
 * Free FFT and convolution (C++)
 *
 * Copyright (c) 2017 Project Nayuki. (MIT License)
 * https://www.nayuki.io/page/free-small-fft-in-multiple-languages
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of
 * this software and associated documentation files (the "Software"), to deal in
 * the Software without restriction, including without limitation the rights to
 * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
 * the Software, and to permit persons to whom the Software is furnished to do so,
 * subject to the following conditions:
 * - The above copyright notice and this permission notice shall be included in
 *   all copies or substantial portions of the Software.
 * - The Software is provided "as is", without warranty of any kind, express or
 *   implied, including but not limited to the warranties of merchantability,
 *   fitness for a particular purpose and noninfringement. In no event shall the
 *   authors or copyright holders be liable for any claim, damages or other
 *   liability, whether in an action of contract, tort or otherwise, arising from,
 *   out of or in connection with the Software or the use or other dealings in the
 *   Software.
 */

#pragma once

#include <complex>
#include <vector>


namespace Fft {
    
    /*
     * Computes the discrete Fourier transform (DFT) of the given complex vector, storing the result back into the vector.
     * The vector can have any length. This is a wrapper function.
     */
    void transform(std::vector<std::complex<double> > &vec);
    
    
    /*
     * Computes the inverse discrete Fourier transform (IDFT) of the given complex vector, storing the result back into the vector.
     * The vector can have any length. This is a wrapper function. This transform does not perform scaling, so the inverse is not a true inverse.
     */
    void inverseTransform(std::vector<std::complex<double> > &vec);
    
    
    /*
     * Computes the discrete Fourier transform (DFT) of the given complex vector, storing the result back into the vector.
     * The vector's length must be a power of 2. Uses the Cooley-Tukey decimation-in-time radix-2 algorithm.
     */
    void transformRadix2(std::vector<std::complex<double> > &vec);
    
    
    /*
     * Computes the discrete Fourier transform (DFT) of the given complex vector, storing the result back into the vector.
     * The vector can have any length. This requires the convolution function, which in turn requires the radix-2 FFT function.
     * Uses Bluestein's chirp z-transform algorithm.
     */
    void transformBluestein(std::vector<std::complex<double> > &vec);
    
    
    /*
     * Computes the circular convolution of the given complex vectors. Each vector's length must be the same.
     */
    void convolve(const std::vector<std::complex<double> > &vecx,const std::vector<std::complex<double> > &vecy,
std::vector<std::complex<double> > &vecout);
    
}

#endif /* SpectrographGen_hpp */

