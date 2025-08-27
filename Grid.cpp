#include "Grid.hpp"


Grid::Grid(int length, int width, int cellSize){
    this->length = length;
    this->width = width;
    this->numRows = width / cellSize;
    this->numColumns = length / cellSize;
    this->numCells = numRows * numColumns;

    for (int i = 0; i < numRows; i++){
        for (int j = 0; j < numColumns; j++){
            gridCells.push_back(Cell(i, i, j, cellSize));
        }
    }
}