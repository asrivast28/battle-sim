#include <iostream>
#include <unistd.h>
#include <tclap/CmdLine.h>

#include "battle_field.hpp"

int main(int argc, char *argv[])
{

    try {
    // define commandline usage
    TCLAP::CmdLine cmd("Collect statistics about a battle.");
    TCLAP::ValueArg<int> iterArg("i", "iterations", "Number of iterations to run", false, 1, "num");
    cmd.add(iterArg);
    TCLAP::ValueArg<int> sizexArg("x", "size-x", "Number of rows.", true, 20, "num");
    cmd.add(sizexArg);
    TCLAP::ValueArg<int> sizeyArg("y", "size-y", "Number of columns.", true, 20, "num");
    cmd.add(sizeyArg);
    TCLAP::ValueArg<int> soldierArg("s", "soldiers", "Rows of soldiers on each side.", false, 3, "num");
    cmd.add(soldierArg);
    TCLAP::ValueArg<int> skill0Arg("g", "skill-green", "Skill of army 0.", false, 150, "num");
    cmd.add(skill0Arg);
    TCLAP::ValueArg<int> skill1Arg("b", "skill-blue", "Skill of army 1.", false, 150, "num");
    cmd.add(skill1Arg);
    cmd.parse(argc, argv);

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

    BattleField::setDynamicFieldDecayFactor(1);
    BattleField::setFollowPreviousProbability(0.25);
    BattleField::setExtendedNeighborhoodSize(5);

    srand(0);
    for (int iter = 0; iter < iterArg.getValue(); ++iter) {


    BattleField ff(sizexArg.getValue(),sizeyArg.getValue());

    // set some soldiers of either type
    for (int i = 0; i < sizeyArg.getValue(); ++i)
    {
        // two rows of soldiers
        /*
        gsoldiers.at(1, i) = Soldier(0, Soldier::SWORDSMAN, randomAtMax(max_skill), randomAtMax(max_aggression));
        gsoldiers.at(2, i) = Soldier(0, Soldier::SWORDSMAN, randomAtMax(max_skill), randomAtMax(max_aggression));
        gsoldiers.at(q*10-2, i) = Soldier(1, Soldier::SWORDSMAN, randomAtMax(max_skill), randomAtMax(max_aggression));
        gsoldiers.at(q*10-1, i) = Soldier(1, Soldier::SWORDSMAN, randomAtMax(max_skill), randomAtMax(max_aggression));
        */
        for (int j = 0; j < soldierArg.getValue(); ++j)
        {
            ff.mat().at(j, i) = Soldier(1, Soldier::SWORDSMAN, skill0Arg.getValue(), 100);
            ff.mat().at(sizexArg.getValue()-j-1, i) = Soldier(0, Soldier::SWORDSMAN, skill1Arg.getValue(), 100);
        }
    }

    ff.setTarget(BattleField::ANNIHILATE_ENEMY);
    // set flags into the middle
    ff.setFlag(0, sizexArg.getValue()/2, sizeyArg.getValue()/2);
    ff.setFlag(1, sizexArg.getValue()/2, sizeyArg.getValue()/2);

    // MUST be called before the start of iterations
    ff.initializeNeighborCounts();

    int i = 0;
    while (true)
    {
        ff.move();
        ff.kill();
        if (ff.status() == BattleField::WON || ff.status() == BattleField::TIED)
        {
            size_t count0 = ff.soldierCount(0);
            size_t count1 = ff.soldierCount(1);
            std::cout << i << "\t" << count0 << "\t" << count1 << std::endl;
            break;
        }
        ++i;
    }
    }
    } catch (TCLAP::ArgException& e) {
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    }
    return 0;
}
