#include "Cell.hpp"


Cell::Cell(int cID, int row, int col, int cellSize)
{
    this->cellID = cID;
    this->cellBlock.setSize(sf::Vector2f(cellSize, cellSize));
    this->cellBlock.setOutlineColor(sf::Color::Black);
    this->cellBlock.setOutlineThickness(2.0f);
    this->cellBlock.setPosition(row * cellSize, col * cellSize);
}
