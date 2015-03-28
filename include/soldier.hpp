/**
 * @file    soldier.hpp
 * @author  Ankit Srivastava <asrivast@gatech.edu>
 * @brief   
 *
 * Copyright (c) TODO
 */

#ifndef SOLDIER_H
#define SOLDIER_H

class Soldier {
public:
    // enumeration type for specifying the type of this soldier
    enum Type {
      LEADER,
      SWORDSMAN,
      ARCHER
    };

    // data type for storing skill level
    typedef unsigned char SkillType;
    // data type for storing aggression level
    typedef unsigned char AggressionType; 

public:
    // Default constructor
    Soldier()
      : m_army(), m_type(),
        m_skill(0), m_aggression(0)
    { }

    // Constructor
    Soldier(bool army, Type type)
        : m_army(army), m_type(type),
          m_skill(0), m_aggression(0)
    { }

    // Copy constructor
    Soldier(const Soldier& soldier)
        : m_army(soldier.m_army), m_type(soldier.m_type),
          m_skill(soldier.m_skill), m_aggression(soldier.m_aggression)
    { }

    // Destructor
    ~Soldier()
    { }

private:
    // the army that this soldier belongs to
    bool m_army;
    // type of this soldier
    Type m_type;
    // skill level of this soldier
    SkillType m_skill;
    // aggression level of this soldier
    AggressionType m_aggression;

}; // class Soldier

#endif // SOLDIER_H
