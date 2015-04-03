
#include <iostream>
#include <floor_field.hpp>
#include <unistd.h>


int main(int argc, char *argv[])
{
    FloorField ff(10,10);

    srand(0);

    // set some soldiers of either type
    for (int i = 0; i < 5; ++i)
    {
        ff.mat().at(1, 2*i + 1) = Soldier(0, Soldier::SWORDSMAN);
        ff.mat().at(8, 2*i) = Soldier(1, Soldier::SWORDSMAN);
    }

    int n_iter = 20;
    for (int i = 0; i < n_iter; ++i)
    {
        ff.move();
        std::cout << ff.mat() << std::endl;
        double sleep_secs = 0.5;
        usleep((unsigned int) (sleep_secs * 1000000));
    }

    return 0;
}
