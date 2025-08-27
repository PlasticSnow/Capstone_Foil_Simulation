#include <SFML/Graphics.hpp>

#pragma once

class Cell{

    public:

        Cell(int cID, int x, int y, int cellSize);

        int getCellID(){return cellID;}
        
        sf::RectangleShape getCellBlock(){return cellBlock;}

    private:
        float leftVelocity = 0;
        float rightVelocity = 0;
        float upwardVelocity = 0;
        float downwardVelocity = 0;
        int cellID;
        sf::RectangleShape cellBlock;
};