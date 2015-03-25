/**
 * @file    mpi_cellular_automata.hpp
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Implements a class for cellular automata, including the data used
 *          to store it.
 *
 * Copyright (c) TODO
 */

#ifndef MPI_CELLULAR_AUTOMATA_H
#define MPI_CELLULAR_AUTOMATA_H

#include <mpi.h>
#include "cellular_automata.hpp"

#define TAG_ROW_COMM 13
#define TAG_COL_COMM 42

class MpiCellularAutomata : public CellularAutomata
{
private:
    // MPI communicators for my column and row
    MPI_Comm row_comm;
    MPI_Comm col_comm;
    int row_id;
    int row_size;
    int col_id;
    int col_size;

    // column wise communication (communicating rows)
    MPI_Request req_send_top;
    MPI_Request req_send_bottom;
    MPI_Request req_recv_top;
    MPI_Request req_recv_bottom;

    // column wise communication (communicating rows)
    MPI_Request req_send_left;
    MPI_Request req_send_right;
    MPI_Request req_recv_left;
    MPI_Request req_recv_right;

    // The MPI data types for rows and columns
    MPI_Datatype row_data_type;
    MPI_Datatype col_data_type;


public:
    MpiCellularAutomata(std::size_t nrows, std::size_t ncols, MPI_Comm row_comm, MPI_Comm col_comm)
        : CellularAutomata(nrows, ncols), row_comm(row_comm), col_comm(col_comm)
    {
        // rank and size of row
        MPI_Comm_rank(row_comm, &col_id);
        MPI_Comm_size(row_comm, &row_size);
        // rank and size of column
        MPI_Comm_rank(col_comm, &row_id);
        MPI_Comm_size(col_comm, &col_size);

        // create the column custom data type (strided)
        MPI_Type_vector(m+2, 1, n+2, MPI_UINT8_T, &col_data_type);
        MPI_Type_commit(&col_data_type);
        // create row data type (congiuous), this is not necessarily needed,
        // but done due to completeness (both cases are handled identical)
        MPI_Type_contiguous(n, MPI_UINT8_T, &row_data_type);
        MPI_Type_commit(&row_data_type);

        // initialize persistent communication
        init_persistent_comm();
    }

    virtual ~MpiCellularAutomata()
    {
        /* free all the MPI_request objects */
        if (row_id > 0)
        {
            // communicate with the top
            MPI_Request_free(&req_send_top);
            MPI_Request_free(&req_recv_top);
        }
        if (row_id < col_size - 1)
        {
            // communicate with the bottom
            MPI_Request_free(&req_send_bottom);
            MPI_Request_free(&req_recv_bottom);
        }

        // row wise communication (communicating columns)
        if (col_id > 0)
        {
            // communicate with the left
            MPI_Request_free(&req_send_left);
            MPI_Request_free(&req_recv_left);
        }
        if (col_id < row_size - 1)
        {
            // communicate with the right
            MPI_Request_free(&req_send_right);
            MPI_Request_free(&req_recv_right);
        }
    }

    /**
     * @brief Generates the next iteration on the parallel computing grid.
     */
    virtual void nextIteration()
    {
        bool needs_col_comm = true;
        // dispatch column communication
        start_comm_colwise();

        // compute core of core
        // for each row i
        for (std::size_t i = 1; i < m-1; ++i)
        {
            // for each column j
            for (std::size_t j = 1; j < n-1; ++j)
            {
                updateCell(i,j);
            }
            // test for column wise communication and dispatch row wise communication
            if ((i%10 == 0) && needs_col_comm && test_comm_colwise()) {
                start_comm_rowwise();
                needs_col_comm = false;
            }
        }

        // enter the block only if row wise communication was not previously started
        if (needs_col_comm) {
            // wait for column wise communication to finish
            wait_comm_colwise();

            // dispatch row wise communication
            start_comm_rowwise();
        }

        // fill the top and bottom rows of the core
        for (std::size_t j = 1; j < n-1; ++j)
        {
            updateCell(0, j);
            updateCell(m-1, j);
        }

        // wait for the row wise communication to finish
        wait_comm_rowwise();

        // fill the outer most columns
        for (std::size_t i = 0; i < m; ++i)
        {
            updateCell(i, 0);
            updateCell(i, n-1);
        }

        // swap buffers
        copyBuffer();
    }

protected:
    /**
     * @brief Initializes the persistent communication channels
     */
    void init_persistent_comm()
    {
        // column wise communication (communicating rows)
        if (row_id > 0)
        {
            // communicate with the top
            MPI_Send_init(data + allindex2offset(1,1), 1, row_data_type,
                    row_id-1, TAG_COL_COMM, col_comm, &req_send_top);
            MPI_Recv_init(data + allindex2offset(0,1), 1, row_data_type,
                    row_id-1, TAG_COL_COMM, col_comm, &req_recv_top);
        }
        if (row_id < col_size - 1)
        {
            // communicate with the bottom
            MPI_Send_init(data + allindex2offset(m,1), 1, row_data_type,
                    row_id+1, TAG_COL_COMM, col_comm, &req_send_bottom);
            MPI_Recv_init(data + allindex2offset(m+1,1), 1, row_data_type,
                    row_id+1, TAG_COL_COMM, col_comm, &req_recv_bottom);
        }

        // row wise communication (communicating columns)
        if (col_id > 0)
        {
            // communicate with the left
            MPI_Send_init(data + allindex2offset(0,1), 1, col_data_type,
                    col_id-1, TAG_ROW_COMM, row_comm, &req_send_left);
            MPI_Recv_init(data + allindex2offset(0,0), 1, col_data_type,
                    col_id-1, TAG_ROW_COMM, row_comm, &req_recv_left);
        }
        if (col_id < row_size - 1)
        {
            // communicate with the right
            MPI_Send_init(data + allindex2offset(0,n), 1, col_data_type,
                    col_id+1, TAG_ROW_COMM, row_comm, &req_send_right);
            MPI_Recv_init(data + allindex2offset(0,n+1), 1, col_data_type,
                    col_id+1, TAG_ROW_COMM, row_comm, &req_recv_right);
        }
    }

    void start_comm_colwise()
    {
        if (row_id > 0)
        {
            // start send and receives
            MPI_Start(&req_send_top);
            MPI_Start(&req_recv_top);
        }
        if (row_id < col_size - 1)
        {
            // start send and receives
            MPI_Start(&req_send_bottom);
            MPI_Start(&req_recv_bottom);
        }
    }

    void start_comm_rowwise()
    {
        if (col_id > 0)
        {
            // start send and receives
            MPI_Start(&req_send_left);
            MPI_Start(&req_recv_left);
        }
        if (col_id < row_size - 1)
        {
            // start send and receives
            MPI_Start(&req_send_right);
            MPI_Start(&req_recv_right);
        }
    }

    void wait_comm_colwise()
    {
        if (row_id > 0)
        {
            MPI_Wait(&req_send_top, MPI_STATUS_IGNORE);
            MPI_Wait(&req_recv_top, MPI_STATUS_IGNORE);
        }
        if (row_id < col_size - 1)
        {
            MPI_Wait(&req_send_bottom, MPI_STATUS_IGNORE);
            MPI_Wait(&req_recv_bottom, MPI_STATUS_IGNORE);
        }
    }

    bool test_comm_colwise()
    {
        int f1=1, f2=1, f3=1, f4=1;
        if (row_id > 0)
        {
            MPI_Test(&req_send_top, &f1, MPI_STATUS_IGNORE);
            MPI_Test(&req_recv_top, &f2, MPI_STATUS_IGNORE);
        }
        if (row_id < col_size - 1)
        {
            MPI_Test(&req_send_bottom, &f3, MPI_STATUS_IGNORE);
            MPI_Test(&req_recv_bottom, &f4, MPI_STATUS_IGNORE);
        }
        return (f1 && f2 && f3 && f4);
    }

    void wait_comm_rowwise()
    {
        if (col_id > 0)
        {
            MPI_Wait(&req_send_left, MPI_STATUS_IGNORE);
            MPI_Wait(&req_recv_left, MPI_STATUS_IGNORE);
        }
        if (col_id < row_size - 1)
        {
            MPI_Wait(&req_send_right, MPI_STATUS_IGNORE);
            MPI_Wait(&req_recv_right, MPI_STATUS_IGNORE);
        }
    }

    bool test_comm_rowwise()
    {
        int f1=1, f2=1, f3=1, f4=1;
        if (col_id > 0)
        {
            MPI_Test(&req_send_left, &f1, MPI_STATUS_IGNORE);
            MPI_Test(&req_recv_left, &f2, MPI_STATUS_IGNORE);
        }
        if (col_id < row_size - 1)
        {
            MPI_Test(&req_send_right, &f3, MPI_STATUS_IGNORE);
            MPI_Test(&req_recv_right, &f4, MPI_STATUS_IGNORE);
        }
        return (f1 && f2 && f3 && f4);
    }
};

#endif // MPI_CELLULAR_AUTOMATA_H
