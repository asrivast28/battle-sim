%module battlesim
%{
#define SWIG_FILE_WITH_INIT
#include "soldier.hpp"
#include "floor_field.hpp"
%}

%include "numpy.i"
%init %{
import_array();
%}

%numpy_typemaps(unsigned char, NPY_UBYTE, std::size_t)
%apply (std::size_t DIM1, std::size_t DIM2, unsigned char* INPLACE_ARRAY2) {(std::size_t, std::size_t, unsigned char*)}

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
  Soldier();
  Soldier(unsigned char, Type, SkillType, AggressionType);
  Soldier(const Soldier&);
};


%nodefaultctor FloorField;
class FloorField {
public:
  FloorField(std::size_t, std::size_t);
  FloorField(std::size_t, std::size_t, unsigned char*);
  void move();
  void kill();
  void setSoldier(std::size_t, std::size_t, const Soldier&);
  void setTarget(const unsigned char, const std::size_t, const std::size_t);
  void initializeNeighborhood();
  void getSoldiers(std::size_t, std::size_t, unsigned char*);
  void printGrid();
  ~FloorField();
};
%clear (std::size_t, std::size_t, unsigned char*);
