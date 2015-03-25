#include <cstdlib>
#include <math.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <mpi.h>

#include "timer.hpp"
#include "mpi_utils.hpp"
#include "mpi_cellular_automata.hpp"


/* Setup timers */
timer counter;
double t_mpi_gol_begin, t_mpi_gol_end;

/**
 * @brief Arranges the processors into a two dimensional grid.
 *
 * The number of processors has to be a perfect square or a power of two.
 *
 * @param p     The number of processors to divide into a grid.
 * @param dims  The output: the number of processors in each dimension.
 */
void get_grid_dimensions(int p, int* dims)
{
    int sqrt_p = (int)sqrt(p);
    // if p is a perfect square root of p
    if (sqrt_p * sqrt_p == p)
    {
        dims[0] = sqrt_p;
        dims[1] = sqrt_p;
        return;
    }
    // if the input is a power of two
    bool powerOfTwo = !(p == 0) && !(p & (p - 1));
    if (powerOfTwo)
    {
        int zero_bits = __builtin_ctz((unsigned int)p);
        int floor_z = zero_bits / 2;
        int ceil_z = zero_bits / 2 + (zero_bits % 2);

        dims[0] = 1 << ceil_z;
        dims[1] = 1 << floor_z;
        return;
    }

    // now we fail
    std::cerr << "Number of processors: " << p << " is not a perfect square of a power of two." << std::endl;
    exit(EXIT_FAILURE);
}

/**
 * @brief Performs a block decomposition.
 *
 * @param n  The number of global elements.
 * @param p  The number of processors.
 * @param i  My local processor rank.
 *
 * @return  The number of local elements.
 */
int get_block_decomposition_n(int n, int p, int i)
{
    int m = n / p;
    int rem = n % p;
    if (i < rem)
        return m + 1;
    else
        return m;
}

/**
 * @brief   Loads the input from file and returns it as big vector.
 *
 * @param filepath              The path to the file to load.
 * @param global_input_params   Additional parameters that are passed out.
 *
 * @return A vector containing the grid from the file.
 */
std::vector<CellularAutomata::cell_t> load_from_file(const char * filepath, int* global_input_params)
{
    std::ifstream infile_stream(filepath, std::ifstream::in);
    std::istream_iterator<int> input_iterator(infile_stream);

    // get number of elements as first element from stream
    int m = *(input_iterator++);
    int n = *(input_iterator++);
    int iterations = *(input_iterator++);

    // assign those to the output params
    global_input_params[0] = m;
    global_input_params[1] = n;
    global_input_params[2] = iterations;

    std::vector<CellularAutomata::cell_t> result(m*n);
    copy_n(input_iterator, m*n, result.begin());

    return result;
}

/**
 * @brief Print the given vector as an output grid of '0' and '1's
 *
 * @param vec    The grid to print as std::vector.
 * @param ncols  The number of columns in the grid.
 * @param os     The output stream to print into (default std::cout)
 */
void printgrid(std::vector<CellularAutomata::cell_t>& vec, int ncols = 0, std::ostream& os = std::cout)
{
    for (unsigned int i = 0; i < vec.size(); ++i)
    {
        if (ncols != 0 && i % ncols == 0)
            os << std::endl;
        os << (int)vec[i] << " ";
    }
    os << std::endl;
}

/**
 * Prints the usage of the program.
 */
void print_usage()
{
    std::cerr << "Usage: ./mpi_gol [options] [input_file]" << std::endl;
    std::cerr << "      Optional arguments:" << std::endl;
    std::cerr << "          -o <file>    Output all solutions to the given file." << std::endl;
    std::cerr << "          -r           Run random grid tests, random grids are generated only locally, no bottleneck at startup." << std::endl;
    std::cerr << "          -x <n>       Sets the number of rows in GoL grid (mandatory with option -r)" << std::endl;
    std::cerr << "          -y <n>       Sets the number of cols in GoL grid (mandatory with option -r)" << std::endl;
    std::cerr << "          -i <n>       Sets the number of iterations for GoL grid (mandatory with option -r)" << std::endl;
    std::cerr << "      Example:" << std::endl;
    std::cerr << "          ./mpi_gol -o output.txt input.txt" << std::endl;
    std::cerr << "                  Plays GoL as specified by input.txt and outputs the result to output.txt" << std::endl;
}

int main(int argc, char *argv[])
{

    /***************************
     *  Parse input arguments  *
     ***************************/

    // forget about first argument (which is the executable's name)
    argc--;
    argv++;

    bool do_output = false;
    char* outfile_path = NULL;

    bool do_read_from_file = true;
    bool do_generate_local = false;
    int gol_row_size = -1, gol_col_size = -1, gol_iterations=-1;
    char* infile_path = NULL;

    // parse optional parameters
    while(argc > 0 && argv[0][0] == '-')
    {
        char option = argv[0][1];
        switch (option)
        {
            case 'o':
                do_output = true;
                // next parameter is file name
                outfile_path = argv[1];
                argv++;
                argc--;
                break;
            case 'r':
                do_generate_local = true;
                do_read_from_file = false;
                break;
            case 'x':
                // the next argument must be the number
                argv++;
                argc--;
                gol_row_size = atoi(argv[0]);
                break;
            case 'y':
                // the next argument must be the number
                argv++;
                argc--;
                gol_col_size = atoi(argv[0]);
                break;
            case 'i':
                // the next argument must be the number
                argv++;
                argc--;
                gol_iterations = atoi(argv[0]);
                break;
            default:
                print_usage();
                exit(EXIT_FAILURE);
        }
        // iterate to next argument
        argv++;
        argc--;
    }

    // check that the mandatory parameters are present
    if (do_read_from_file && argc < 1)
    {
        print_usage();
        exit(EXIT_FAILURE);
    }
    if (do_generate_local && (gol_row_size<0 || gol_col_size<0 || gol_iterations<0))
    {
        print_usage();
        exit(EXIT_FAILURE);
    }

    // parse mandatory parameters
    if (do_read_from_file)
    {
        infile_path = argv[0];
    }

    /********************
     *  Setting up MPI  *
     ********************/
    // set up MPI
    int p;
    // int rank;
    MPI_Init(&argc, &argv);
    // get total size of processors and current rank
    //MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    if (p > 1) {
        // get the dimensions of the processor grid
        int dims[2];
        get_grid_dimensions(p, dims);

        int periods[2] = {0, 0};
        MPI_Comm grid_comm;
        MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &grid_comm);
        // now get column and row communicators
        MPI_Comm row_comm;
        MPI_Comm col_comm;
        int row_dims[2] = {0, 1};
        int col_dims[2] = {1, 0};
        MPI_Cart_sub(grid_comm, row_dims, &row_comm);
        MPI_Cart_sub(grid_comm, col_dims, &col_comm);

        // get our column index and row index
        int col_id;
        MPI_Comm_rank(row_comm, &col_id);
        int row_id;
        MPI_Comm_rank(col_comm, &row_id);

        int global_input_params[3] = {gol_row_size, gol_col_size, gol_iterations};
        std::vector<CellularAutomata::cell_t> row_input;
        if (do_read_from_file) {
            // distribute all the data from rank 0 processor to first column
            if (col_id == 0)
            {
                std::vector<CellularAutomata::cell_t> global_input;
                if (row_id == 0)
                {
                    // read the file into global array
                    global_input = load_from_file(infile_path, global_input_params);
                }

                // broadcast global parameters across first column
                MPI_Bcast(global_input_params, 3, MPI_INT, 0, col_comm);

                // distribute data across first column
                row_input = scatter_vector_block_decomp(global_input, col_comm);
            }
        }

        if ((do_generate_local) && (col_id == 0)) {
            // broadcast global parameters across first column
            MPI_Bcast(global_input_params, 3, MPI_INT, 0, col_comm);
        }
        // broadcast global parameters across the rows
        MPI_Bcast(global_input_params, 3, MPI_INT, 0, row_comm);

        int num_iterations = global_input_params[2];

        // distribute data across all rows
        std::vector<CellularAutomata::cell_t> local_input;
        if (do_read_from_file) {
            local_input = scatter_rows_block_decomp(row_input, global_input_params[1], row_comm);
        }

        int local_nrows = get_block_decomposition_n(global_input_params[0], dims[0], row_id);
        int local_ncols = get_block_decomposition_n(global_input_params[1], dims[1], col_id);

        //std::cout << "processor = (" << row_id << ", " << col_id << std::endl;
        //std::cout << "Rows and cols: " << local_nrows << " " << local_ncols << " " << num_iterations << std::endl;
        //printgrid(local_input, local_ncols);

        // create and initialize the local automata
        MpiCellularAutomata localAutomata(local_nrows, local_ncols, row_comm, col_comm);
        if (do_generate_local) {
            localAutomata.initRandom();
        } else {
            localAutomata.init(local_input);
        }

        // run the iterations
        // Timing only on master node (and barrierized)
        MPI_Barrier (grid_comm);
        t_mpi_gol_begin = counter.get_ms();
        for (int iter = 0; iter < num_iterations; ++iter)
        {
            localAutomata.nextIteration();
        }
        MPI_Barrier (grid_comm);
        t_mpi_gol_end = counter.get_ms();

        if (do_output) {
            // gather all the data again
            local_input = localAutomata.getCells();

            // gather rows into first column
            row_input = gather_rows(local_input, local_nrows, row_comm);

            // within the first column, gather everything to processor (0,0)
            if (col_id == 0)
            {
                std::vector<CellularAutomata::cell_t> global_input;
                global_input = gather_vectors(row_input, col_comm);

                // if i am processor (0,0), then output result
                if (row_id == 0)
                {
                    if (outfile_path != NULL)
                    {
                        // print to file
                        std::ofstream outstr(outfile_path);
                        printgrid(global_input, global_input_params[1], outstr);
                        outstr.close();
                    }
                    //else
                    //{
                        //printgrid(global_input, global_input_params[1]);
                    //}
                }
            }
        }
    } else {
        //Run the sequential algorithm on one processor
        int global_input_params[3] = {gol_row_size, gol_col_size, gol_iterations};
        std::vector<CellularAutomata::cell_t> global_input;
        if (do_read_from_file) {
            global_input = load_from_file(infile_path, global_input_params);
        }

        // create and initialize the local automata
        CellularAutomata localAutomata(global_input_params[0], global_input_params[1]);
        if (do_generate_local) {
            localAutomata.initRandom();
        } else {
            localAutomata.init(global_input);
        }

        // run the iterations
        // Timing
        t_mpi_gol_begin = counter.get_ms();
        for (int iter = 0; iter < global_input_params[2]; ++iter)
        {
            localAutomata.nextIteration();
        }
        t_mpi_gol_end = counter.get_ms();

        if (do_output) {
            // gather all the data again
            global_input = localAutomata.getCells();
            if (outfile_path != NULL)
            {
                // print to file
                std::ofstream outstr(outfile_path);
                printgrid(global_input, global_input_params[1], outstr);
                outstr.close();
            }
            //else
            //{
                //printgrid(global_input, global_input_params[1]);
            //}
        }
    }

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        std::cerr << "gol took: " << (t_mpi_gol_end - t_mpi_gol_begin) << " ms" << std::endl;
    }

    // clean up mpi
    MPI_Finalize();

    return 0;
}
