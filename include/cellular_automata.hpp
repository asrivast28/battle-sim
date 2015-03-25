/**
 * @file    cellular_automata.hpp
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Implements a class for cellular automata, including the data used
 *          to store it.
 *
 * Copyright (c) TODO
 */

#ifndef CELLULAR_AUTOMATA_H
#define CELLULAR_AUTOMATA_H


#include <cstdlib>
#include <assert.h>
#include <algorithm>


class CellularAutomata
{
public:
    /// define the base data type for a cell
    typedef uint8_t cell_t;
protected:
    // the allocated data for double buffering
    cell_t * data;
    cell_t * back_buffer;

    /// Number of rows
    const std::size_t m;
    /// Number of columns
    const std::size_t n;
    /// Total data size
    const std::size_t data_size;

    // returns the offset for the given 2d indeces relative to the core
    inline std::size_t index2offset(const std::size_t i, const std::size_t j) const
    {
        assert(0 <= i && i < m);
        assert(0 <= j && j < n);
        return (i+1) * (n+2) + (j+1);
    }

    // returns the offset for the given 2d indeces relative to the whole data
    inline std::size_t allindex2offset(const std::size_t i, const std::size_t j) const
    {
        assert(0 <= i && i < m+2);
        assert(0 <= j && j < n+2);
        return i * (n+2) + j;
    }

    /// returns the value of the cell given by the indeces relative to the core
    inline cell_t getCell(const std::size_t i, const std::size_t j) const
    {
        return *(data + index2offset(i,j));
    }

    /// returns the value of the cell given by the indeces
    inline cell_t getCellGlobal(const std::size_t i, const std::size_t j) const
    {
        return *(data + allindex2offset(i,j));
    }

    /// sets the value for the given cell indeces in the back-buffer.
    inline void setBufferCell(const std::size_t i, const std::size_t j, const cell_t value)
    {
        *(back_buffer + index2offset(i,j)) = value;
    }

    ///'sets the value for the given cell indeces in the front-buffer
    inline void setCell(const std::size_t i, const std::size_t j, const cell_t value)
    {
        *(data + index2offset(i,j)) = value;
    }

    /// copies the back buffer into the front buffer
    inline void copyBuffer()
    {
        std::copy(back_buffer, back_buffer+data_size, data);
    }

    /// Game of Life
    inline void updateCell(std::size_t i, std::size_t j)
    {
        // TODO: this could potentially be done more efficiently
        // count number of neighbors with status = 1
        int count = 0;

        // convert to global index:
        i++;
        j++;

        // NW
        count += getCellGlobal(i-1,j-1);
        // N
        count += getCellGlobal(i-1,j);
        // NE
        count += getCellGlobal(i-1,j+1);

        // W
        count += getCellGlobal(i, j-1);
        // E
        count += getCellGlobal(i, j+1);

        // SW
        count += getCellGlobal(i+1,j-1);
        // S
        count += getCellGlobal(i+1,j);
        // SE
        count += getCellGlobal(i+1,j+1);

        // if this cell is empty
        i--;
        j--;
        if (getCell(i,j) == 0)
        {
            if (count == 3)
            {
                // become occupied
                setBufferCell(i,j,1);
            }
            else
                setBufferCell(i,j,0);
        } else {
            if (count == 2 || count == 3)
            {
                // stay occupied
                setBufferCell(i,j,1);
            }
            else
            {
                setBufferCell(i,j,0);
            }
        }
    }

public:

    /// Construtor taking two arguments for the size of the grid
    CellularAutomata(std::size_t nrows, std::size_t ncols)
        : m(nrows), n(ncols), data_size((m+2)*(n+2))
    {
        // allocate the data buffers
        data = new cell_t[data_size];
        back_buffer = new cell_t[data_size];

        // initialize everything with zero
        std::fill(data, data+data_size, 0);
        std::fill(back_buffer, back_buffer+data_size, 0);
    }

    /// Deconstructor
    virtual ~CellularAutomata()
    {
        // deallocate the data buffers
        delete[] data;
        data = NULL;
        delete[] back_buffer;
        back_buffer = NULL;
    }

    /**
     * @brief Sequentially calculates the next iteration of the cellular automata.
     */
    virtual void nextIteration()
    {
        // for each cell:

        // for each row i
        for (std::size_t i = 0; i < m; ++i)
        {
            // for each column j
            for (std::size_t j = 0; j < n; ++j)
            {
                updateCell(i,j);
            }
        }
        // swap buffers
        copyBuffer();
    }

    /**
     * @brief   Initializes the cellular automata grid with the given values.
     *
     * @param input_data
     */
    void init(std::vector<cell_t>& input_data)
    {
        // initializes the inner grid with the given data
        assert(m*n == input_data.size());

        // fill in row major order
        std::vector<cell_t>::iterator input_it = input_data.begin();
        for (unsigned int i = 0; i < m; ++i)
        {
            for (unsigned int j = 0; j < n; ++j)
            {
                /* code */
                setCell(i, j, *(input_it++));
            }
        }
    }

    /**
     * @brief Randomly initializes the 
     */
    void initRandom()
    {
        for (unsigned int i = 0; i < m; ++i)
        {
            for (unsigned int j = 0; j < n; ++j)
            {
                /* code */
                setCell(i, j, rand() % 2);
            }
        }
    }

    /**
     * @brief Returns a vector with values of all cells.
     *
     * @return 
     */
    std::vector<cell_t> getCells()
    {
        std::vector<cell_t> result(m*n);
        // fill in row major order
        std::vector<cell_t>::iterator out_it = result.begin();
        for (unsigned int i = 0; i < m; ++i)
        {
            for (unsigned int j = 0; j < n; ++j)
            {
                /* code */
                *(out_it++) = getCell(i, j);
            }
        }
        return result;
    }

    inline std::size_t getColumnStride()
    {
        return n+2;
    }

    inline std::size_t getNumColumns()
    {
        return n;
    }

    inline std::size_t getNumRows()
    {
        return m;
    }
};


#endif // CELLULAR_AUTOMATA_H
