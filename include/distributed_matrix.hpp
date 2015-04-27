/**
 * @file    distributed_matrix.hpp
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Templated matrix and border_matrix types.
 *
 * Copyright (c) 2014 Georgia Institute of Technology. All Rights Reserved.
 */
#ifndef DISTRIBUTED_MATRIX_HPP
#define DISTRIBUTED_MATRIX_HPP

#include <mpi.h>
#include <assert.h>

#include "matrix.hpp"
#include "datatypes.hpp"
#include <prettyprint.hpp>

#define TAG_COL_COMM 13
#define TAG_ROW_COMM 42

template <typename T>
class distributed_matrix : public border_matrix<T>
{
public:
    // assume local sizes are all equal
    matrix<T> gather() {
        // gather the inner parts of the border matrizes to (0,0)
        matrix<T> mat;
        int rank, p;
        MPI_Comm_rank(comm, &rank);
        MPI_Comm_size(comm, &p);
        mxx::datatype<T> dt;
        MPI_Datatype send_type, rcv_type;
        MPI_Type_vector(this->nrows(), this->ncols(), this->ncols()+2*this->border_size, dt.type(), &send_type);
        MPI_Type_commit(&send_type);
        MPI_Type_contiguous(this->nrows()*this->ncols(), dt.type(), &rcv_type);
        MPI_Type_commit(&rcv_type);

        // 1) gather to first column
        std::vector<T> rowdata;
        if (rank_y == 0) {
            rowdata.resize(size_y*this->ncols()*this->nrows());
            MPI_Gather(&this->at(0,0), 1, send_type, &rowdata[0], 1, rcv_type, 0, row_comm);
        } else {
            MPI_Gather(&this->at(0,0), 1, send_type, NULL, 0, MPI_DATATYPE_NULL, 0, row_comm);
        }

        // 2) gather to (0,0)
        if (rank_y == 0) {
            // local reorder
            std::vector<T> buf(rowdata.size());
            for (int i = 0; i < size_y; ++i) {
                for (std::size_t x = 0; x < this->nrows(); x++) {
                    for (std::size_t y = 0; y < this->ncols(); y++) {
                        buf[x*(this->ncols()*size_y) + i*this->ncols() + y] = rowdata[i*(this->nrows()*this->ncols()) + x*this->ncols() + y];
                    }
                }
            }
            // then gather to (0,0)
            int size = this->nrows()*this->ncols()*size_y;
            if (rank_x == 0) {
                mat = matrix<T>(size_x*this->nrows(), size_y*this->ncols());
                MPI_Gather(&buf[0], size, dt.type(), &mat.at(0,0), size, dt.type(), 0, col_comm);
            } else {
                MPI_Gather(&buf[0], size, dt.type(), NULL, 0, MPI_DATATYPE_NULL, 0, col_comm);
            }
        }
        return mat;
    }

    void scatter(matrix<T> mat) {
        int rank, p;
        MPI_Comm_rank(comm, &rank);
        MPI_Comm_size(comm, &p);
        mxx::datatype<T> dt;
        MPI_Datatype send_type, rcv_type;
        MPI_Type_vector(this->nrows(), this->ncols(), this->nrows()+2*this->border_size, dt.type(), &send_type);
        MPI_Type_commit(&send_type);
        MPI_Type_contiguous(this->nrows()*this->ncols(), dt.type(), &rcv_type);
        MPI_Type_commit(&rcv_type);

        // 2) scatter from (0,0)
        std::vector<T> rowdata;
        if (rank_y == 0) {
            // local reorder
            rowdata.resize(size_y*this->ncols()*this->nrows());
            std::vector<T> buf(rowdata.size());
            // then gather to (0,0)
            int size = this->nrows()*this->ncols()*size_y;
            if (rank_x == 0) {
                MPI_Scatter(&mat.at(0,0), size, dt.type(), &buf[0], size, dt.type(), 0, col_comm);
            } else {
                MPI_Scatter(NULL, 0, MPI_DATATYPE_NULL, &buf[0], size, dt.type(), 0, col_comm);
            }
            for (int i = 0; i < size_y; ++i) {
                for (std::size_t x = 0; x < this->nrows(); x++) {
                    for (std::size_t y = 0; y < this->ncols(); y++) {
                        rowdata[i*(this->nrows()*this->ncols()) + x*this->ncols() + y] = buf[x*(this->ncols()*size_y) + i*this->ncols() + y];
                    }
                }
            }
        }
        // 1) gather to first column
        if (rank_y == 0) {
            //rowdata.resize(size_y*this->ncols()*this->nrows());
            MPI_Scatter(&rowdata[0], 1, rcv_type, &this->at(0,0), 1, send_type, 0, row_comm);
        } else {
            MPI_Scatter(NULL, 0, MPI_DATATYPE_NULL, &this->at(0,0), 1, send_type, 0, row_comm);
        }

    }

    void sync_boundaries() {
        // send inner elements and receive into boundaries

        // create the column custom data type (strided)
        MPI_Datatype col_data_type;
        mxx::datatype<T> dt;
        MPI_Type_vector(this->nrows(), this->border_size, this->ncols()+2*this->border_size, dt.type(), &col_data_type);
        MPI_Type_commit(&col_data_type);
        MPI_Datatype row_dt;
        MPI_Type_contiguous((this->ncols()+2*this->border_size)*this->border_size, dt.type(), &row_dt);
        MPI_Type_commit(&row_dt);

        std::vector<MPI_Request> reqs;
        reqs.reserve(4);
        // 1) across rows send columns
        if (rank_y > 0) {
            // recv from left
            reqs.push_back(MPI_Request());
            MPI_Irecv(&this->at(0,-this->border_size), 1, col_data_type, rank_y-1, TAG_ROW_COMM, row_comm, &reqs.back());
            // send to left
            reqs.push_back(MPI_Request());
            MPI_Isend(&this->at(0, 0), 1, col_data_type, rank_y-1, TAG_ROW_COMM, row_comm, &reqs.back());
        }
        if (rank_y < size_y-1) {
            // recv from right
            reqs.push_back(MPI_Request());
            MPI_Irecv(&this->at(0,this->ncols()), 1, col_data_type, rank_y+1, TAG_ROW_COMM, row_comm, &reqs.back());
            // send to right
            reqs.push_back(MPI_Request());
            MPI_Isend(&this->at(0, this->ncols()-this->border_size), 1, col_data_type, rank_y+1, TAG_ROW_COMM, row_comm, &reqs.back());
        }
        MPI_Waitall(reqs.size(), &reqs[0], MPI_STATUSES_IGNORE);
        reqs.clear(); reqs.reserve(4);


        // 2) across columns sending rows
        if (rank_x > 0) {
            // recv from top
            reqs.push_back(MPI_Request());
            MPI_Irecv(&this->at(-this->border_size,-this->border_size), 1, row_dt, rank_x-1, TAG_COL_COMM, col_comm, &reqs.back());
            // send to top
            reqs.push_back(MPI_Request());
            MPI_Isend(&this->at(0,-this->border_size), 1, row_dt, rank_x-1, TAG_COL_COMM, col_comm, &reqs.back());
        }
        if (rank_x < size_x-1) {
            // recv from bottom
            reqs.push_back(MPI_Request());
            MPI_Irecv(&this->at(this->nrows(),-this->border_size), 1, row_dt, rank_x+1, TAG_COL_COMM, col_comm, &reqs.back());
            // send to bottom
            reqs.push_back(MPI_Request());
            MPI_Isend(&this->at(this->nrows()-this->border_size, -this->border_size), 1, row_dt, rank_x+1, TAG_COL_COMM, col_comm, &reqs.back());
        }
        MPI_Waitall(reqs.size(), &reqs[0], MPI_STATUSES_IGNORE);
        reqs.clear(); reqs.reserve(4);

    }


    template <typename L>
    void accumulate_back(L func) {
        // communicate back (from boundaries into the regular cell and
        // accumulate there)
        MPI_Datatype col_data_type;
        mxx::datatype<T> dt;
        MPI_Type_vector(this->nrows(), this->border_size, this->ncols()+2*this->border_size, dt.type(), &col_data_type);
        MPI_Type_commit(&col_data_type);

        MPI_Datatype col_rcv_type;
        MPI_Type_contiguous(this->nrows()*this->border_size, dt.type(), &col_rcv_type);
        MPI_Type_commit(&col_rcv_type);

        MPI_Datatype row_dt;
        MPI_Type_contiguous((this->ncols()+2*this->border_size)*this->border_size, dt.type(), &row_dt);
        MPI_Type_commit(&row_dt);

        std::vector<MPI_Request> reqs;
        reqs.reserve(4);

        // 1) across columns sending rows
        std::vector<T> buf1((this->ncols()+2*this->border_size)*this->border_size);
        std::vector<T> buf2((this->ncols()+2*this->border_size)*this->border_size);
        if (rank_x > 0) {
            // recv from top
            reqs.push_back(MPI_Request());
            MPI_Irecv(&buf1[0], 1, row_dt, rank_x-1, TAG_COL_COMM, col_comm, &reqs.back());
            // send to top
            reqs.push_back(MPI_Request());
            MPI_Isend(&this->at(-this->border_size,-this->border_size), 1, row_dt, rank_x-1, TAG_COL_COMM, col_comm, &reqs.back());
        }
        if (rank_x < size_x-1) {
            // recv from bottom
            reqs.push_back(MPI_Request());
            MPI_Irecv(&buf2[0], 1, row_dt, rank_x+1, TAG_COL_COMM, col_comm, &reqs.back());
            // send to bottom
            reqs.push_back(MPI_Request());
            MPI_Isend(&this->at(this->nrows(),-this->border_size), 1, row_dt, rank_x+1, TAG_COL_COMM, col_comm, &reqs.back());
        }
        MPI_Waitall(reqs.size(), &reqs[0], MPI_STATUSES_IGNORE);
        reqs.clear(); reqs.reserve(4);

        // accumulate buffers into active area
        if (rank_x > 0) {
            index_t i = 0;
            for (index_t x = 0; x < this->border_size; ++x) {
                for (index_t y = -this->border_size; y < (index_t)(this->ncols() + this->border_size); ++y) {
                    this->at(x, y) = func(this->at(x, y), buf1[i++]);
                }
            }
        }
        if (rank_x < size_x-1) {
            index_t i = 0;
            for (index_t x = this->nrows() - this->border_size; x < this->nrows(); ++x) {
                for (index_t y = -this->border_size; y < (index_t)(this->ncols() + this->border_size); ++y) {
                    this->at(x, y) = func(this->at(x, y), buf2[i++]);
                }
            }
        }


        buf1.resize(this->nrows()*this->border_size);
        buf2.resize(this->nrows()*this->border_size);
        // 1) across rows send columns
        if (rank_y > 0) {
            // recv from left
            reqs.push_back(MPI_Request());
            MPI_Irecv(&buf1[0], 1, col_rcv_type, rank_y-1, TAG_ROW_COMM, row_comm, &reqs.back());
            // send to left
            reqs.push_back(MPI_Request());
            MPI_Isend(&this->at(0,-this->border_size), 1, col_data_type, rank_y-1, TAG_ROW_COMM, row_comm, &reqs.back());
        }
        if (rank_y < size_y-1) {
            // recv from right
            reqs.push_back(MPI_Request());
            MPI_Irecv(&buf2[0], 1, col_rcv_type, rank_y+1, TAG_ROW_COMM, row_comm, &reqs.back());
            // send to right
            reqs.push_back(MPI_Request());
            MPI_Isend(&this->at(0,this->ncols()), 1, col_data_type, rank_y+1, TAG_ROW_COMM, row_comm, &reqs.back());
        }
        MPI_Waitall(reqs.size(), &reqs[0], MPI_STATUSES_IGNORE);

        // accumulate buffers into active area
        if (rank_y > 0) {
            index_t i = 0;
            for (index_t x = 0; x < this->nrows(); ++x) {
                for (index_t y = 0; y < this->border_size; ++y) {
                    this->at(x, y) = func(this->at(x, y), buf1[i++]);
                }
            }
        }
        if (rank_y < size_y-1) {
            index_t i = 0;
            for (index_t x = 0; x < this->nrows(); ++x) {
                for (index_t y = (index_t)this->ncols()-this->border_size; y < (index_t)this->ncols(); ++y) {
                    this->at(x, y) = func(this->at(x, y), buf2[i++]);
                }
            }
        }

        // and sync changes back
        sync_boundaries();
    }

    distributed_matrix() : border_matrix<T>() {}

    distributed_matrix(const distributed_matrix<T>& o)
        : border_matrix<T>(o),
          comm(o.comm), row_comm(o.row_comm), col_comm(o.col_comm),
          size_x(o.size_x), size_y(o.size_y),
          rank_x(o.rank_x), rank_y(o.rank_y) {}

    distributed_matrix(std::size_t nrows, std::size_t ncols, std::size_t bordersize, MPI_Comm comm)
        : border_matrix<T>(nrows, ncols, bordersize), comm(comm) {
        // assume comm is a 2D communicator
        int dims[2];    // sizes of each dimension
        int periods[2]; // perioids of each dimension (should be {0,0})
        int coords[2];  // coordinates of this processor in the grid
        MPI_Cart_get(comm, 2, dims, periods, coords);
        size_x = dims[0];
        size_y = dims[1];
        rank_x = coords[0];
        rank_y = coords[1];
        int row_dims[2] = {0, 1};
        int col_dims[2] = {1, 0};
        MPI_Cart_sub(comm, row_dims, &row_comm);
        MPI_Cart_sub(comm, col_dims, &col_comm);
    }

    distributed_matrix(std::size_t nrows, std::size_t ncols, std::size_t bordersize, MPI_Comm comm, const T& init)
        : border_matrix<T>(nrows, ncols, bordersize, init), comm(comm) {
        // assume comm is a 2D communicator
        int dims[2];    // sizes of each dimension
        int periods[2]; // perioids of each dimension (should be {0,0})
        int coords[2];  // coordinates of this processor in the grid
        MPI_Cart_get(comm, 2, dims, periods, coords);
        size_x = dims[0];
        size_y = dims[1];
        rank_x = coords[0];
        rank_y = coords[1];
        int row_dims[2] = {0, 1};
        int col_dims[2] = {1, 0};
        MPI_Cart_sub(comm, row_dims, &row_comm);
        MPI_Cart_sub(comm, col_dims, &col_comm);
    }

    virtual ~distributed_matrix() {}

private:
    // cartesian 2d grid
    MPI_Comm comm;
    MPI_Comm row_comm;
    MPI_Comm col_comm;
    // coordinates of processor
    int rank_x;
    int rank_y;
    int size_x;
    int size_y;
};



#endif // MATRIXH_HPP
