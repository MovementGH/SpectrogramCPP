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
#include <thread>
class SpectrographGen {
public:
    SpectrographGen();
    SpectrographGen(int SpecSampleRate,int SpecResPerSample,float Compression);
    bool loadFromFile(std::string File);
    void saveToFile(std::string File);
protected:
    sf::Image generateSpectogram();
    void generateLine(sf::Image* Output,int LineStart,int LineEnd);
    void sort(std::vector<std::complex<float>>& tmp,std::complex<float>* a,int n);
    void fft(std::vector<std::complex<float>>& tmp,std::complex<float>* x,int n,int s);
    std::vector<sf::Int16> m_Samples;
    sf::Int32 m_SampleRate;
    int m_SpecSampleRate;
    int m_SpecResPerSample;
    float m_Compression;
    int m_FFTSize;
    std::vector<float> m_Hanning;
    std::vector<std::vector<std::complex<float>>> m_Polar;
};

class SpectrographDecode {
public:
    SpectrographDecode();
    SpectrographDecode(int SpecSampleRate,int SampleRate,int ChannelCount,float Compression);
    bool loadFromFile(std::string File);
    void saveToFile(std::string File);
protected:
    sf::SoundBuffer generateBuffer();
    void decodeLine(sf::Int16* Output,int StartLine,int EndLine);
    void sort(std::vector<std::complex<float>>& tmp,std::complex<float>* a,int n);
    void fft(std::vector<std::complex<float>>& tmp,std::complex<float>* x,int n,int s);
    void ifft(std::vector<std::complex<float>>& tmp,std::vector<std::complex<float>>& x,int s);
    sf::Image m_Image;
    sf::Int32 m_SampleRate;
    int m_SpecSampleRate;
    int m_SpecResPerSample;
    int m_ChannelCount;
    float m_Compression;
    int m_FFTSize;
    std::vector<float> m_Hanning;
    std::vector<std::vector<std::complex<float>>> m_Polar;
};
#endif
