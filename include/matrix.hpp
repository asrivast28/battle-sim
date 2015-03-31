/**
 * @file    matrix.hpp
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Templated matrix and border_matrix types.
 *
 * Copyright (c) 2014 Georgia Institute of Technology. All Rights Reserved.
 */
#ifndef MATRIXH_HPP
#define MATRIXH_HPP

#include <assert.h>
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
        return *this;
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

    /// element access operator with `at(i,j)` (otherwise identical)
    inline const T& at(std::size_t i, std::size_t j) const
    {
        assert(0 <= i && i < m);
        assert(0 <= j && j < n);
        return m_data[i*n + j];
    }

    /// element access operator with `at(i,j)` (otherwise identical)
    inline T& at(std::size_t i, std::size_t j)
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
    /// number of rows
    std::size_t m;
    /// number of columns
    std::size_t n;
    /// data of the matrix
    std::vector<T> m_data;
};

// matrix with border (generalized)
// TODO the floor field model will consist of multiple of these!
template <typename T>
class border_matrix : protected matrix<T> {
public:
    typedef T value_type;
    typedef matrix<T> matrix_type;

    /// Default constructor
    border_matrix() : matrix_type(), border_size(0) {}

    /// Copy constructor
    border_matrix(const border_matrix<T>& o)
        : matrix_type(o),
          border_size(o.border_size) {}

    /// emtpy initialization
    border_matrix(std::size_t nrows, std::size_t ncols, std::size_t border_size)
        : matrix_type(nrows+2*border_size, ncols+2*border_size), border_size(border_size) {}

    /// initialization with data from iterators
    template <typename InputIterator>
    border_matrix(std::size_t nrows, std::size_t ncols, std::size_t border_size, InputIterator begin, InputIterator end)
        : matrix_type(nrows+2*border_size, ncols+2*border_size), border_size(border_size)
    {
        std::size_t size = std::distance(begin, end);
        if (size == nrows*ncols)
        {
            // only data for internal matrix given -> have to copy row by row
            for (std::size_t i = 0; i < nrows; ++i)
            {
                std::copy(begin+i*ncols, begin + (i+1)*ncols,
                          &this->m_data[i*(this->n) + border_size]);
            }
        }
        else if (size == this->n * this->m)
        {
            // data for the whole matrix, including border is given
            std::copy(begin, end, this->m_data.begin());
        }
        else
        {
            assert(false);
        }
    }

    /// initialization with data
    border_matrix(std::size_t nrows, std::size_t ncols, std::size_t border_size, const std::vector<T>& data)
        : border_matrix(nrows, ncols, border_size, data.begin(), data.end()) {
    }

    /// Destructor
    virtual ~border_matrix() {}

    /// Assignment operator
    virtual border_matrix<T>& operator=(const border_matrix<T>& other) {
        this->m = other.m;
        this->n = other.n;
        this->m_data = other.m_data;
        border_size = other.border_size;
        return *this;
    }

    /// element access operator
    inline const T& operator()(std::size_t i, std::size_t j) const
    {
        assert(0 <= i && i < this->m-2*border_size);
        assert(0 <= j && j < this->n-2*border_size);
        return this->m_data[(i+border_size)*this->n + j + border_size];
    }

    /// element access operator
    inline T& operator()(std::size_t i, std::size_t j)
    {
        assert(0 <= i && i < this->m-2*border_size);
        assert(0 <= j && j < this->n-2*border_size);
        return this->m_data[(i+border_size)*this->n + j + border_size];
    }

    /// element access operator with `at(i,j)` (otherwise identical)
    inline const T& at(std::size_t i, std::size_t j) const
    {
        assert(0 <= i && i < this->m-2*border_size);
        assert(0 <= j && j < this->n-2*border_size);
        return this->m_data[(i+border_size)*this->n + j + border_size];
    }

    /// element access operator with `at(i,j)` (otherwise identical)
    inline T& at(std::size_t i, std::size_t j)
    {
        assert(0 <= i && i < this->m-2*border_size);
        assert(0 <= j && j < this->n-2*border_size);
        return this->m_data[(i+border_size)*this->n + j + border_size];
    }

    /// returns a raw pointer to the underlying data
    inline T* data() {
        return &this->m_data[0];
    }

    /// Returns the matrix including border
    inline matrix_type& mat() {
        return *((matrix_type*) this);
    }

    /// Returns the number of rows
    inline std::size_t nrows() const {
        return this->m - 2*border_size;
    }

    /// Returns the number of columns
    inline std::size_t ncols() const {
        return this->n - 2*border_size;
    }

    /// Returns the size of the border
    inline std::size_t bordersize() const {
        return border_size;
    }

protected:
    std::size_t border_size;
};


#endif // MATRIXH_HPP
