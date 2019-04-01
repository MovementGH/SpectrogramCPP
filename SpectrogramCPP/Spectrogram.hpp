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
    void sort(std::complex<float>* tmp,std::complex<float>* a,int n);
    void fft(std::complex<float>* tmp,std::complex<float>* x,int n,int s);
    void ifft(std::complex<float>* tmp,std::vector<std::complex<float>>& x,int s);
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
    void generateLine(sf::Uint8* Output,int LineStart,int LineEnd,sf::Vector2i& Size);
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
