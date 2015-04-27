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
#include <iostream>

// using signed index type (for border matrix)
typedef int64_t index_t;

// simple matrix type, data saved as row-major
template <typename T>
class matrix {
public:
    typedef T value_type;

    /// Default constructor
    matrix() : m(0), n(0), m_data() {}

    /// Copy constructor
    matrix(const matrix<T>& o) : m(o.m), n(o.n), m_data(o.m_data) {}

    /// empty initialization
    matrix(std::size_t nrows, std::size_t ncols)
        : m(nrows), n(ncols), m_data(n*m) {}

    /// value initialization
    matrix(std::size_t nrows, std::size_t ncols, const T& init)
        : m(nrows), n(ncols), m_data(n*m, init) {}

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
    inline const T& operator()(index_t i, index_t j) const
    {
        assert(0 <= i && i < m);
        assert(0 <= j && j < n);
        return m_data[i*n + j];
    }

    /// element access operator
    inline T& operator()(index_t i, index_t j)
    {
        assert(0 <= i && i < m);
        assert(0 <= j && j < n);
        return m_data[i*n + j];
    }

    /// element access operator with `at(i,j)` (otherwise identical)
    inline const T& at(index_t i, index_t j) const
    {
        assert(0 <= i && i < m);
        assert(0 <= j && j < n);
        return m_data[i*n + j];
    }

    /// element access operator with `at(i,j)` (otherwise identical)
    inline T& at(index_t i, index_t j)
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
    inline std::size_t  ncols() const {
        return n;
    }

    /// set everything to zero or default constructed
    inline void zero() {
        for (index_t i = 0; i < m; ++i)
            for (index_t j = 0; j < n; ++j)
                at(i,j) = T();
    }
protected:
    /// number of rows
    index_t m;
    /// number of columns
    index_t n;
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

    border_matrix(std::size_t nrows, std::size_t ncols, std::size_t border_size, const T& init)
        : matrix_type(nrows+2*border_size, ncols+2*border_size, init), border_size(border_size) {}

    /// initialization with data from iterators
    template <typename InputIterator>
    border_matrix(std::size_t nrows, std::size_t ncols, std::size_t border_size, InputIterator begin, InputIterator end)
        : matrix_type(nrows+2*border_size, ncols+2*border_size), border_size(border_size)
    {
        std::size_t size = std::distance(begin, end);
        if (size == nrows*ncols)
        {
            // only data for internal matrix given -> have to copy row by row
            for (index_t i = 0; i < this->m; ++i)
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
    border_matrix(index_t nrows, index_t ncols, index_t border_size, const std::vector<T>& data)
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
    inline const T& operator()(index_t i, index_t j) const
    {
        assert(-border_size <= i && i < this->m-border_size);
        assert(-border_size <= j && j < this->n-border_size);
        return this->m_data[(i+border_size)*this->n + j + border_size];
    }

    /// element access operator
    inline T& operator()(index_t i, index_t j)
    {
        assert(-border_size <= i && i < this->m-border_size);
        assert(-border_size <= j && j < this->n-border_size);
        return this->m_data[(i+border_size)*this->n + j + border_size];
    }

    /// element access operator with `at(i,j)` (otherwise identical)
    inline const T& at(index_t i, index_t j) const
    {
        assert(-border_size <= i && i < this->m-border_size);
        assert(-border_size <= j && j < this->n-border_size);
        return this->m_data[(i+border_size)*this->n + j + border_size];
    }

    /// element access operator with `at(i,j)` (otherwise identical)
    inline T& at(index_t i, index_t j)
    {
        assert(-border_size <= i && i < this->m-border_size);
        assert(-border_size <= j && j < this->n-border_size);
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

    /// set everything to zero or default constructed
    inline void zero() {
        this->mat().zero();
    }
    // set the borders to zero
    inline void border_zero() {
        for (index_t i = 0; i < this->m; ++i)
            for (index_t j = 0; j < this->n; ++j)
                if ((i < border_size || i >= (index_t)this->nrows()+border_size) ||
                    (j < border_size || j >= (index_t)this->ncols()+border_size))
                    mat().at(i,j) = T();
    }

protected:
    index_t border_size;
};


template<typename T>
std::ostream& operator<< (std::ostream& stream, const matrix<T>& mat)
{
    for (std::size_t i = 0; i < mat.nrows(); ++i)
    {
        stream << "[";
        for (std::size_t j = 0; j < mat.ncols(); ++j)
        {
            if (j != 0)
                stream << " ";
            stream << mat(i,j);
        }
        stream << "]" << std::endl;
    }
    return stream;
}

template<typename T>
std::ostream& operator<< (std::ostream& stream, const border_matrix<T>& mat)
{
    for (std::size_t i = 0; i < mat.nrows(); ++i)
    {
        stream << "[";
        for (std::size_t j = 0; j < mat.ncols(); ++j)
        {
            if (j != 0)
                stream << " ";
            stream << mat(i,j);
        }
        stream << "]" << std::endl;
    }
    return stream;
}

template<>
std::ostream& operator<< (std::ostream& stream, const matrix<unsigned char>& mat)
{
    for (std::size_t i = 0; i < mat.nrows(); ++i)
    {
        stream << "[";
        for (std::size_t j = 0; j < mat.ncols(); ++j)
        {
            if (j != 0)
                stream << " ";
            stream << static_cast<unsigned>(mat(i,j));
        }
        stream << "]" << std::endl;
    }
    return stream;
}

#endif // MATRIXH_HPP
