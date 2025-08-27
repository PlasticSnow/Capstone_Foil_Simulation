#include <SFML/Graphics.hpp>
#include <iostream>
#include <string>
#include "Grid.hpp"

int main() {
    int winSizeX = 500;
    int winSizeY = 500;

    Grid* grid = new Grid(winSizeX, winSizeY, 25);

    

    sf::RenderWindow window(sf::VideoMode(winSizeX, winSizeY), "SFML Test");

    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();
        }

        window.clear(sf::Color::Cyan);

        // Draw Grids Squares
        for (int i = 0; i < grid->getNumGridCells(); i++){
            window.draw(grid->getGridCell(i).getCellBlock());
        }

        window.display();
    }

    return 0;
}
