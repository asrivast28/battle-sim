#include <mpi.h>
#include <cmath>
#include <cstdlib>
#include <unistd.h>

#include <distributed_matrix.hpp>

int main(int argc, char *argv[])
{
    // set up MPI
    MPI_Init(&argc, &argv);

    // get communicator size and my rank
    MPI_Comm comm = MPI_COMM_WORLD;
    int p, rank;
    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &rank);

    /* code */
    /* ... */

    int q = (int)sqrt(p);
    if (p != q*q) {
        throw std::runtime_error("The number of MPI processes must be a perfect square");
    }
    MPI_Comm grid_comm;
    int dims[2] = {q, q};
    int periods[2] = {0, 0};
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &grid_comm);

    // blah
    distributed_matrix<int> dmat(4, 4, 2, grid_comm);

    matrix<int> mat;
    if (rank == 0)
    {
        mat = matrix<int>(q*4, q*4);
        for (unsigned int i = 0; i < mat.nrows(); ++i)
            for (unsigned int j = 0; j < mat.ncols(); ++j)
                mat(i,j) = rand() % 2;
        std::cout << "Master mat:" << std::endl;
        std::cout << mat << std::endl;
    }
    dmat.scatter(mat);
    std::cout << "mat on p=" << rank << std::endl;
    std::cout << dmat.mat() << std::endl;
    dmat.sync_boundaries();
    //usleep(rank*1000000);
    //std::cout << dmat << std::endl;
    std::cout << dmat.mat() << std::endl;
    dmat.accumulate_back(std::plus<int>());
    std::cout  << dmat.mat() << std::endl;

    // finalize MPI
    MPI_Finalize();
    return 0;
}
