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
class FFTBase {
public:
    FFTBase(int FFTWidth);
    void init(int FFTWidth);
protected:
    void sort(std::vector<std::complex<float>>& tmp,std::complex<float>* a,int n);
    void fft(std::vector<std::complex<float>>& tmp,std::complex<float>* x,int n,int s);
    void ifft(std::vector<std::complex<float>>& tmp,std::vector<std::complex<float>>& x,int s);
    std::vector<float> m_Hanning;
    std::vector<std::vector<std::complex<float>>> m_Polar;
    int m_FFTSize;
    int m_FFTWidth;
};
class SpectrographGen:public FFTBase {
public:
    SpectrographGen();
    SpectrographGen(int SpecSampleRate,int SpecResPerSample,float Compression);
    bool loadFromFile(std::string File);
    void loadFromBuffer(sf::SoundBuffer& Buffer);
    void saveToFile(std::string File);
    sf::Image generateImage();
protected:
    void generateLine(sf::Image* Output,int LineStart,int LineEnd);
    std::vector<sf::Int16> m_Samples;
    sf::Int32 m_SampleRate;
    int m_SpecSampleRate;
    float m_Compression;
};
class SpectrographDecode:public FFTBase {
public:
    SpectrographDecode();
    SpectrographDecode(int SpecSampleRate,int SampleRate,int ChannelCount,float Compression);
    bool loadFromFile(std::string File);
    void loadFromImage(sf::Image& Image);
    void saveToFile(std::string File);
    sf::SoundBuffer generateBuffer();
protected:
    void decodeLine(sf::Int16* Output,int StartLine,int EndLine);
    sf::Image m_Image;
    sf::Int32 m_SampleRate;
    int m_SpecSampleRate;
    int m_ChannelCount;
    float m_Compression;
};
#endif
