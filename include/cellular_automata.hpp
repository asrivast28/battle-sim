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
#include <vector>


// simple matrix type, data saved as row-major
template <typename T>
class matrix {
public:
    typedef T value_type;

    /// Default constructor
    matrix() : m(0), n(0), m_data() {}

    /// Copy constructor
    matrix(const matrix<T>& o) : m(o.m), n(o.n), m_data(o.m_data) {}

    /// emtpy initialization
    matrix(std::size_t nrows, std::size_t ncols)
        : m(nrows), n(ncols), m_data(n*m) {}

    /// initialization with data
    matrix(std::size_t nrows, std::size_t ncols, const std::vector<T>& data)
        : m(nrows), n(ncols), m_data(data) {}

    /// initialization with data from iterators
    template <typename InputIterator>
    matrix(std::size_t nrows, std::size_t ncols, InputIterator begin, InputIterator end)
        : m(nrows), n(ncols), m_data(begin, end) {
        assert(std::distance(begin, end) == n*m);
    }

    /// Destructor
    virtual ~matrix() {}

    /// Assignment operator
    virtual matrix<T>& operator=(const matrix<T>& other) {
        m = other.m;
        n = other.n;
        m_data = other.m_data;
    }

    /// element access operator
    inline const T& operator()(std::size_t i, std::size_t j) const
    {
        assert(0 <= i && i < m);
        assert(0 <= j && j < n);
        return m_data[i*n + j];
    }

    /// element access operator
    inline T& operator()(std::size_t i, std::size_t j)
    {
        assert(0 <= i && i < m);
        assert(0 <= j && j < n);
        return m_data[i*n + j];
    }

    /// returns a raw pointer to the underlying data
    inline T* data() {
        return &m_data[0];
    }

    /// Returns the number of rows
    inline std::size_t nrows() const {
        return m;
    }

    /// Returns the number of columns
    inline std::size_t ncols() const {
        return n;
    }

protected:
    /// numer rows
    const std::size_t m;
    /// number columns
    const std::size_t n;
    /// data of the matrix
    std::vector<T> m_data;
};

// matrix with border (generalized)
// TODO the floor field model will consist of multiple of these!
// simple matrix type, data saved as row-major
template <typename T>
class border_matrix : protected matrix<T> {
public:
    typedef T value_type;
    typedef matrix<T> matrix_type;

    /// Default constructor
    border_matrix() : matrix_type() {}

    /// Copy constructor
    border_matrix(const border_matrix<T>& o)
        : matrix_type(o),
          border_size(o.border_size) {}

    /// emtpy initialization
    border_matrix(std::size_t nrows, std::size_t ncols, std::size_t border_size)
        : matrix_type(nrows+2*border_size, ncols+2*border_size) {}

    /// TODO!!!
    /// initialization with data from iterators
    template <typename InputIterator>
    border_matrix(std::size_t nrows, std::size_t ncols, InputIterator begin, InputIterator end)
        : m(nrows), n(ncols), m_data(begin, end) {
        assert(std::distance(begin, end) == n*m);
    }

    /// initialization with data
    border_matrix(std::size_t nrows, std::size_t ncols, std::size_t border_size, const std::vector<T>& data)
        : matrix_type(nrows+2*border_size, ncols+2*border_size) {
        // two cases: either the vector gives data for the whole grid,
        // including the border, or only for the internal part
        if (data.size() == nrows*ncols)
        {
            // only data for internal matrix given -> have to copy row by row
            for (std::size_t i = 0; i < nrows; ++i)
            {
                std::copy(&data[i*ncols], &data[(i+1)*ncols],
                          &this->m_data[i*(this->n) + border_size]);
            }
        }
        else if (data.size() == this->n * this->m)
        {
            // data for the whole matrix, including border is given
            this->m_data = data;
        }
    }


    /// Destructor
    virtual ~border_matrix() {}

    /// Assignment operator
    virtual border_matrix<T>& operator=(const border_matrix<T>& other) {
        this->m = other.m;
        this->n = other.n;
        this->m_data = other.m_data;
        border_size = other.border_size;
    }

    /// element access operator
    inline const T& operator()(std::size_t i, std::size_t j) const
    {
        assert(0 <= i && i < this->m);
        assert(0 <= j && j < this->n);
        return m_data[i*n + j];
    }

    /// element access operator
    inline T& operator()(std::size_t i, std::size_t j)
    {
        assert(0 <= i && i < m);
        assert(0 <= j && j < n);
        return m_data[i*n + j];
    }

    /// returns a raw pointer to the underlying data
    inline T* data() {
        return &m_data[0];
    }

    /// Returns the number of rows
    inline std::size_t nrows() const {
        return m;
    }

    /// Returns the number of columns
    inline std::size_t ncols() const {
        return n;
    }

protected:
    std::size_t border_size;
};


template <typename T>
class CellularAutomata {
public:
    /// define the base data type for a cell
    typedef T cell_type;
protected:
    // the allocated data for double buffering
    std::vector<cell_type> data;
    std::vector<cell_type> back_buffer;
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
    inline cell_type getCell(const std::size_t i, const std::size_t j) const
    {
        return data[index2offset(i,j)];
    }

    /// returns the value of the cell given by the indeces
    inline cell_type getCellGlobal(const std::size_t i, const std::size_t j) const
    {
        return data[allindex2offset(i,j)];
    }

    /// sets the value for the given cell indeces in the back-buffer.
    inline void setBufferCell(const std::size_t i, const std::size_t j, const cell_type value)
    {
        back_buffer[index2offset(i,j)] = value;
    }

    ///'sets the value for the given cell indeces in the front-buffer
    inline void setCell(const std::size_t i, const std::size_t j, const cell_type value)
    {
        data[index2offset(i,j)] = value;
    }

    /// copies the back buffer into the front buffer
    inline void copyBuffer()
    {
        std::copy(back_buffer, back_buffer+data_size, data);
    }

    /// Game of Life
    inline void updateCell(std::size_t i, std::size_t j)
    {
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
        : m(nrows), n(ncols), data_size((m+2)*(n+2)),
          data(data_size), back_buffer(data_size)
    {
        // allocate the data buffers
        data = new cell_type[data_size];
        back_buffer = new cell_type[data_size];
    }

    /// Deconstructor
    virtual ~CellularAutomata() {}

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
    void init(std::vector<cell_type>& input_data)
    {
        // initializes the inner grid with the given data
        assert(m*n == input_data.size());

        // fill in row major order
        typename std::vector<cell_type>::iterator input_it = input_data.begin();
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
    std::vector<cell_type> getCells()
    {
        std::vector<cell_type> result(m*n);
        // fill in row major order
        typename std::vector<cell_type>::iterator out_it = result.begin();
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
