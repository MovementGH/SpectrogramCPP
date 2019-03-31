#include <iostream>
#include <SpectrographGen.hpp>
#include <SFML/Audio.hpp>
#define OUTPUTPATH std::string("/users/mayo/desktop/")
int main(int argc, const char * argv[]) {
    sf::Clock Timer;
    SpectrographGen Generator(250,4096,2);
    Generator.loadFromFile("come alive.ogg");
    sf::Image Gram=Generator.generateImage();
    SpectrographDecode Decoder(250,44100,2,2);
    Decoder.loadFromImage(Gram);
    Decoder.saveToFile(OUTPUTPATH+"test.ogg");
    std::cout<<Timer.getElapsedTime().asSeconds();
    return 0;
}
