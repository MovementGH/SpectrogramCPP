//
//  SpectrographGen.cpp
//  LibSpectrographCPP
//
//  Created by Mayo Furgerson on 3/27/19.
//  Copyright © 2019 WimMa Games. All rights reserved.
//

#include "SpectrographGen.hpp"


sf::Image SpectrographGen::generateSpectogram() {
    sf::Image Result;
    Result.create(m_SpecSampleRate*(m_Samples.size()/m_SampleRate),m_SpecResPerSample/8);
    std::vector<std::complex<float>> WorkingSamples(m_SpecResPerSample);
    for(int i=0;i<Result.getSize().x;i++) {
        std::size_t start=std::min((std::size_t)i*(m_SampleRate/m_SpecSampleRate),m_Samples.size()-m_SpecResPerSample-1);
        for(std::size_t i2=start;i2<start+m_SpecResPerSample;i2++)
            WorkingSamples[i2-start]=std::complex<float>(m_Samples[i2]*m_Hanning[i2-start],0);
        fft(WorkingSamples.data(),(int)WorkingSamples.size(),m_FFTSize);
        sf::Uint8 res=0;
        for(int i2=WorkingSamples.size()*.375;i2<WorkingSamples.size()/2;i2++) {
            res=std::min(std::max(WorkingSamples[i2+WorkingSamples.size()/2].real()/85000+128,0.f),255.f);
            if(res>126&&res<130)
                res=128;
            Result.setPixel(i,i2-WorkingSamples.size()*.375,{res,res,res,255});
        }
    }
    return Result;
}
void SpectrographGen::saveToFile(std::string File) {
    generateSpectogram().saveToFile(File);
}
bool SpectrographGen::loadFromFile(std::string File) {
    sf::SoundBuffer Buffer;
    if(Buffer.loadFromFile(File)) {
        m_SampleRate=Buffer.getSampleRate();
        m_Samples.resize(Buffer.getSampleCount());
        for(int i=0;i<Buffer.getSampleCount();i++)
            m_Samples[i]=Buffer.getSamples()[i];
        return true;
    }
    return false;
}
SpectrographGen::SpectrographGen() : m_Samples(), m_SampleRate(0), m_SpecSampleRate(0), m_SpecResPerSample(0) {
    
}
SpectrographGen::SpectrographGen(int SpecSampleRate, int SpecResPerSample) : m_Samples(), m_SampleRate(0), m_SpecSampleRate(SpecSampleRate), m_SpecResPerSample(SpecResPerSample) {
    m_Temp.resize(m_SpecResPerSample/2);
    m_Hanning.resize(m_SpecResPerSample);
    for(int i(0);i<m_SpecResPerSample;i++)
        m_Hanning[i]=0.54-0.46*cos(2*M_PI*i/(float)m_SpecResPerSample);
    m_FFTSize=log2(m_SpecResPerSample);
    m_Polar.resize(m_FFTSize,std::vector<std::complex<float>>(m_SpecResPerSample/2));
    for(int i=0;i<m_FFTSize;i++)
        for(int i2=0;i2<m_SpecResPerSample/2;i2++)
            m_Polar[i][i2]=(exp(std::complex<float>(0,-2*M_PI/(pow(2,i+1))*i2)));
}
void SpectrographGen::sort(std::complex<float>* a,int n) {
    for(int i=0;i<n;i++)
        m_Temp[i]=a[i*2+1];
    for(int i=0;i<n;i++)
        a[i]=a[i*2];
    for(int i=0;i<n;i++)
        a[i+n]=m_Temp[i];
}
void SpectrographGen::fft(std::complex<float>* x,int n,int s) {
    if(n>1) {
        n/=2,sort(&x[0],n),fft(&x[0],n,s-1),fft(&x[0]+n,n,s-1);
        for(int k=0;k<n;k++)
            t=m_Polar[s-1][k]*x[k+n],
            x[k+n]=x[k]-t,
            x[k]+=t;
    }
}
sf::SoundBuffer SpectrographDecode::generateBuffer() {
    std::vector<sf::Int16> Samples(((float)m_Image.getSize().x/(float)m_SpecSampleRate)*m_SampleRate+m_SpecResPerSample);
    std::vector<std::complex<float>> WorkingSamples(m_SpecResPerSample,std::complex<float>(0,0));
    std::vector<sf::Int16> WorkingSamples2(m_SpecResPerSample,0);
    for(int x=0;x<m_Image.getSize().x;x++) {
        for(int i=0;i<WorkingSamples.size()/8;i++)
            WorkingSamples2[WorkingSamples.size()/2+i]=m_Image.getPixel(x,m_Image.getSize().y-i-1).g-128,
            WorkingSamples2[WorkingSamples.size()/2-i]=WorkingSamples2[WorkingSamples.size()/2+i];
        for(int i2=0;i2<WorkingSamples.size();i2++)
            WorkingSamples[i2]={WorkingSamples2[i2]*85000.f,0};
        ifft(WorkingSamples,m_FFTSize);
        for(int i2=WorkingSamples.size()*.03f;i2<WorkingSamples.size()*.97f;i2++)
            Samples[x*(m_SampleRate/m_SpecSampleRate)+i2]=(WorkingSamples[i2].real()/500)*m_Hanning[i2];
    }
    sf::SoundBuffer Buffer;
    Buffer.loadFromSamples(Samples.data(),Samples.size(),m_ChannelCount,m_SampleRate);
    return Buffer;
}
void SpectrographDecode::saveToFile(std::string File) {
    generateBuffer().saveToFile(File);
}
bool SpectrographDecode::loadFromFile(std::string File) {
    if(m_Image.loadFromFile(File)) {
        m_SpecResPerSample=m_Image.getSize().y*8;
        m_Temp.resize(m_SpecResPerSample/2);
        m_Hanning.resize(m_SpecResPerSample);
        for(int i(0);i<m_SpecResPerSample;i++)
            m_Hanning[i]=0.54-0.46*cos(2*M_PI*i/(float)m_SpecResPerSample);
        m_FFTSize=log2(m_SpecResPerSample);
        m_Polar.resize(m_FFTSize,std::vector<std::complex<float>>(m_SpecResPerSample/2));
        for(int i=0;i<m_FFTSize;i++)
            for(int i2=0;i2<m_SpecResPerSample/2;i2++)
                m_Polar[i][i2]=(exp(std::complex<float>(0,-2*M_PI/(pow(2,i+1))*i2)));
        return true;
    }
    return false;
}
SpectrographDecode::SpectrographDecode() : m_SampleRate(0), m_SpecSampleRate(0) {
    
}
SpectrographDecode::SpectrographDecode(int SpecSampleRate, int SampleRate,int ChannelCount) : m_SpecSampleRate(SpecSampleRate), m_SampleRate(SampleRate), m_ChannelCount(ChannelCount) {
    
}
void SpectrographDecode::sort(std::complex<float>* a,int n) {
    for(int i=0;i<n;i++)
        m_Temp[i]=a[i*2+1];
    for(int i=0;i<n;i++)
        a[i]=a[i*2];
    for(int i=0;i<n;i++)
        a[i+n]=m_Temp[i];
}
void SpectrographDecode::fft(std::complex<float>* x,int n,int s) {
    if(n>1) {
        n/=2,sort(&x[0],n),fft(&x[0],n,s-1),fft(&x[0]+n,n,s-1);
        for(int k=0;k<n;k++)
            t=m_Polar[s-1][k]*x[k+n],x[k+n]=x[k]-t,x[k]+=t;
    }
}
void SpectrographDecode::ifft(std::vector<std::complex<float>>& x,int s) {
    std::transform(x.cbegin(),x.cend(),x.begin(),static_cast<std::complex<double>(*)(const std::complex<double>&)>(std::conj));
    fft(x.data(),x.size(),s);
    std::transform(x.cbegin(),x.cend(),x.begin(),static_cast<std::complex<double>(*)(const std::complex<double>&)>(std::conj));
}
