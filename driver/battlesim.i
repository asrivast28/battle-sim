%module battlesim
%{
#define SWIG_FILE_WITH_INIT
#include "soldier.hpp"
#include "battle_field.hpp"
%}

%include <std_map.i>
%include <std_pair.i>
%include <std_vector.i>

%include "numpy.i"
%init %{
import_array();
%}

%numpy_typemaps(unsigned char, NPY_UBYTE, size_t)
%apply (size_t DIM1, size_t DIM2, unsigned char* INPLACE_ARRAY2) {(size_t, size_t, unsigned char*)}

%numpy_typemaps(bool, NPY_BOOL, size_t)
%apply (size_t DIM1, size_t DIM2, bool* INPLACE_ARRAY2) {(size_t, size_t, bool*)}

void seed(int);

class Soldier {
public:
    enum Type {
        LEADER,
        SWORDSMAN,
        ARCHER,
        EMPTY
    };

    typedef unsigned char SkillType;
    typedef unsigned char AggressionType;

public:
    static void setKillRadiusMap(const std::map<Soldier::Type, unsigned char>&);
    static void setDynamicFieldMap(const std::map<Soldier::Type, unsigned char>&);

public:
    Soldier();
    Soldier(const unsigned char, const Type, const SkillType, const AggressionType);
    Soldier(const Soldier&);
};

%template() std::pair<Soldier::Type, unsigned char>;
%template(SoldierTypeMap) std::map<Soldier::Type, unsigned char>;

%template() std::pair<size_t, unsigned char>;
%template() std::pair<size_t, Soldier>;
%template(KilledPositionVector) std::vector<size_t>;
%template(SoldierPositionVector) std::vector<std::pair<size_t, unsigned char> >;
%template(SoldierVector) std::vector<std::pair<size_t, Soldier> >;

%nodefaultctor BattleField;
class BattleField {
public:
    enum Target {
        CAPTURE_FLAG,
        ANNIHILATE_ENEMY
    };

    enum Status {
        ONGOING,
        WON,
        TIED
    };

public:
    static void setExtendedNeighborhoodSize(const unsigned char);
    static void setDynamicFieldDecayFactor(const unsigned char);
    static void setFollowPreviousProbability(const float);

public:
    BattleField(const size_t, const size_t);
    BattleField(size_t, size_t, unsigned char*);
    void setTarget(const Target);
    void setSoldiers(const std::vector<std::pair<size_t, Soldier> >&);
    void setFlag(const unsigned char, const size_t, const size_t);
    void move(std::vector<std::pair<size_t, unsigned char> >&);
    size_t kill(std::vector<size_t>&, bool = true);
    size_t kill();
    size_t soldierCount(const unsigned char) const;
    Status status() const;
    unsigned char winner() const;
    void printGrid();
    ~BattleField();
};
%clear (size_t, size_t, unsigned char*);
