#include <vector>
#include "Cell.hpp"

#pragma include

class Grid{

    public:

        Grid(int length, int width, int cellSize);

        Cell getGridCell(int cellIndex){return gridCells[cellIndex];}

        int getNumGridCells(){return numCells;}

    private:
        int length;
        int width;
        int numCells;
        int numRows;
        int numColumns;
        std::vector<Cell> gridCells;;

};