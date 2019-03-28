//
//  main.cpp
//  SpectrographCPP
//
//  Created by Mayo Furgerson on 3/27/19.
//  Copyright © 2019 WimMa Games. All rights reserved.
//

#include <iostream>
#include <SpectrographGen.hpp>
#include <SFML/Audio.hpp>

int main(int argc, const char * argv[]) {
    // insert code here...
    sf::Clock Timer;
//    SpectrographGen Generator(200,4096);
//    Generator.loadFromFile("come alive.ogg");
//    Generator.saveToFile("/users/mayo/desktop/test.png");
    SpectrographDecode Decoder(200,44100,2);
    Decoder.loadFromFile("/users/mayo/desktop/test.png");
    Decoder.saveToFile("/users/mayo/desktop/test.ogg");
    std::cout<<Timer.getElapsedTime().asSeconds();
    return 0;
}
