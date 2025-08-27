#include <SFML/Graphics.hpp>

#pragma once

class Cell{

    public:

        Cell(int cID, int row, int col, int cellSize);

        int getCellID(){return cellID;}
        
        sf::RectangleShape getCellBlock(){return cellBlock;}


        /* Cell Edge Velocity Getters */ 
        float getLeftVel(){return leftVelocity;}
        
        float getRightVel(){return rightVelocity;}
        
        float getTopVel(){return topVelocity;}
        
        float getBottomVel(){return bottomVelocity;}

        /* Cell Edge Velocity Setters */
        void setLeftVel(float velocity){leftVelocity = velocity;}

        void setRightVel(float velocity){rightVelocity = velocity;}

        void setTopVel(float velocity){topVelocity = velocity;}

        void setBottomVel(float velocity){bottomVelocity = velocity;}


    private:

        // Velocity of Cell's Left Edge
        float leftVelocity = 0;

        // Velocity of Cell's Right Edge
        float rightVelocity = 0;

        // Velocity of Cell's Top Edge
        float topVelocity = 0;

        // Velocity of Cell's Bottom Edge
        float bottomVelocity = 0;

        int cellID;

        sf::RectangleShape cellBlock;
};