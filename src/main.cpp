
#include <iostream>
#include <unistd.h>

#include "battle_field.hpp"

int main(int argc, char *argv[])
{
    BattleField ff(20,20);

    srand(0);

    Soldier::SkillType max_skill = 0;
    max_skill = ~max_skill;
    Soldier::AggressionType max_aggression = 0;
    max_aggression = ~max_aggression;

    // set some soldiers of either type
    for (int i = 0; i < 20; ++i)
    {
        ff.mat().at(0, i) = Soldier(0, Soldier::SWORDSMAN, randomAtMax(max_skill), randomAtMax(max_aggression));
        ff.mat().at(1, i) = Soldier(0, Soldier::SWORDSMAN, randomAtMax(max_skill), randomAtMax(max_aggression));
        ff.mat().at(18, i) = Soldier(1, Soldier::SWORDSMAN, randomAtMax(max_skill), randomAtMax(max_aggression));
        ff.mat().at(19, i) = Soldier(1, Soldier::SWORDSMAN, randomAtMax(max_skill), randomAtMax(max_aggression));
    }
    ff.setFlag(0, 0, 5);
    ff.setFlag(1, 10, 5); 

    // MUST be called before the start of iterations
    //ff.initializeNeighborCounts();

    int n_iter = 10;
    double sleep_secs = 0.5;
    std::cout << ff.mat() << std::endl;
    usleep((unsigned int) (sleep_secs * 1000000));
    for (int i = 0; i < n_iter; ++i)
    {
        std::cout << "------------------------------------------" << std::endl;

        ff.move();
        std::cout << ff.mat() << std::endl;
        usleep((unsigned int) (sleep_secs * 1000000));

        ff.kill();
        std::cout << ff.mat() << std::endl;
        usleep((unsigned int) (sleep_secs * 1000000));

    }

    return 0;
}
