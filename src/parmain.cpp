
#include <mpi.h>
#include <iostream>
#include <unistd.h>

#include "par_battle_field.hpp"

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


    srand(13317);

    if (argc < 2) {
        std::cout << "Usage: ./parbattle <global_size>" << std::endl;
    }

    int size = atoi(argv[1]);
    //int sizex = 200;
    //int sizey = 200;
    int sizex = size/q;
    int sizey = size/q;
    int rows = 10;

    // configure
    Soldier::SkillType max_skill = 0;
    max_skill = ~max_skill;
    Soldier::AggressionType max_aggression = 0;
    max_aggression = ~max_aggression;
    std::map<Soldier::Type, unsigned char> killRadiusMap;
    killRadiusMap[Soldier::LEADER] = 1;
    killRadiusMap[Soldier::SWORDSMAN] = 1;
    killRadiusMap[Soldier::ARCHER] = 5;
    Soldier::setKillRadiusMap(killRadiusMap);

    ParallelBattleField::setDynamicFieldDecayFactor(1);
    ParallelBattleField::setFollowPreviousProbability(0.25);
    ParallelBattleField::setExtendedNeighborhoodSize(4);

    ParallelBattleField ff(sizex, sizey, grid_comm);

    matrix<Soldier> gsoldiers;
    if (rank == 0) {
        // set some soldiers of either type
        gsoldiers = matrix<Soldier>(q*sizex, q*sizey);
        for (int i = 0; i < q*sizey; ++i)
        {
            // two rows of soldiers
            /*
            gsoldiers.at(1, i) = Soldier(0, Soldier::SWORDSMAN, randomAtMax(max_skill), randomAtMax(max_aggression));
            gsoldiers.at(2, i) = Soldier(0, Soldier::SWORDSMAN, randomAtMax(max_skill), randomAtMax(max_aggression));
            gsoldiers.at(q*10-2, i) = Soldier(1, Soldier::SWORDSMAN, randomAtMax(max_skill), randomAtMax(max_aggression));
            gsoldiers.at(q*10-1, i) = Soldier(1, Soldier::SWORDSMAN, randomAtMax(max_skill), randomAtMax(max_aggression));
            */
            for (int j = 0; j < rows; ++j)
            {
                gsoldiers.at(j, i) = Soldier(0, Soldier::SWORDSMAN, 255, 100);
                gsoldiers.at(q*sizex-j-1, i) = Soldier(1, Soldier::SWORDSMAN, 100, 100);
            }
        }
    }
    ff.mat().scatter(gsoldiers);
    //ff.setFlag(0, q*sizex-1, sizey*q/2);
    //ff.setFlag(1, 0, sizey*q/2);
    ff.setTarget(ParallelBattleField::ANNIHILATE_ENEMY);
    ff.setFlag(0, sizex*q/2, sizey*q/2);
    ff.setFlag(1, sizex*q/2, sizey*q/2);

    // MUST be called before the start of iterations
    //ff.initializeNeighborCounts();

    //std::cout << ff.mat() << std::endl;
    matrix<Soldier> mat = ff.mat().gather();
    //std::cout << mat << std::endl;


    //int n_iter = 30;
    double sleep_secs = 0.2;
    //usleep((unsigned int) (sleep_secs * 1000000));
    //for (int i = 0; i < n_iter; ++i)
    int i = 0;
    while (true)
    {
        //std::cout << "------------------------------------------" << std::endl;
        ff.move();
        //if (rank == 0)
        size_t count0 = ff.globalSoldierCount(0);
        size_t count1 = ff.globalSoldierCount(1);
        if (i % 10 == 0 && rank == 0)
            std::cout << "at iteration i = " << i << " surivors: " << count0 << ", " << count1 << std::endl;
        //if (rank == 0)
        //    std::cout << "Army 0: " << count0 << ", Army 1: " << count1 << std::endl;
        //matrix<Soldier> gmat = ff.mat().gather();
        //if (rank == 0)
        //    std::cout << gmat << std::endl;
        //std::cout << ff.mat() << std::endl;
        //usleep((unsigned int) (sleep_secs * 1000000));

        ff.kill();
        //gmat = ff.mat().gather();
        //if (rank == 0)
        //    std::cout << gmat << std::endl;
        //std::cout << ff.mat() << std::endl;
        //usleep((unsigned int) (sleep_secs * 1000000));

        if (ff.status() == ParallelBattleField::WON || ff.status() == ParallelBattleField::TIED)
        {
            if (ff.status() == ParallelBattleField::TIED) {
                if (rank == 0)
                    std::cout << "Tied!!!!" << std::endl;
            }
            else if (rank == 0)
                    std::cout << "Winner: " << (unsigned)ff.winner() << " at iteration i = " << i << std::endl;
            break;
        }
        ++i;
    }

    MPI_Finalize();
    return 0;
}
