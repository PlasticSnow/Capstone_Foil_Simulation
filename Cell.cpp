#include "Cell.hpp"


Cell::Cell(int cID, int x, int y, int cellSize)
{
    this->cellID = cID;
    this->cellBlock.setSize(sf::Vector2f(cellSize, cellSize));
    this->cellBlock.setOutlineColor(sf::Color::Black);
    this->cellBlock.setOutlineThickness(2.0f);
    this->cellBlock.setPosition(x * cellSize, y * cellSize);
}
