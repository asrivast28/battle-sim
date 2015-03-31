/**
 * @file    soldier.hpp
 * @author  Ankit Srivastava <asrivast@gatech.edu>
 * @brief
 *
 * Copyright (c) TODO
 */

#ifndef SOLDIER_H
#define SOLDIER_H

#include <map>
#include <iostream>

class Soldier {
public:
    // enumeration type for specifying the type of this soldier
    enum Type {
      LEADER,
      SWORDSMAN,
      ARCHER,
      EMPTY
    };

    // data type for storing skill level
    typedef unsigned char SkillType;
    // data type for storing aggression level
    typedef unsigned char AggressionType;

public:
    /// default constructor
    Soldier()
      : m_army(), m_type(EMPTY),
        m_skill(0), m_aggression(0)
    { }

    /// constructor
    Soldier(bool army, Type type)
        : m_army(army), m_type(type),
          m_skill(0), m_aggression(0)
    { }

    /// copy constructor
    Soldier(const Soldier& soldier)
        : m_army(soldier.m_army), m_type(soldier.m_type),
          m_skill(soldier.m_skill), m_aggression(soldier.m_aggression)
    { }

    /// returns index of the army to which this soldier belongs
    unsigned char
    army() const
    {
        return m_army ? 1 : 0;
    }

    /// returns index of the enemy army
    unsigned char
    enemy() const
    {
        return !m_army ? 1 : 0;
    }

    // returns if this cell is empty
    bool empty() const
    {
        return m_type == EMPTY;
    }

    void kill()
    {
        m_type = EMPTY;
    }

    /// returns aggression of this soldier as a fraction of the maximum aggression
    float
    aggression() const
    {
        AggressionType max_aggresion = 0;
        max_aggresion = ~max_aggresion;
        return static_cast<float>(m_aggression) / max_aggresion;
    }

    /// returns skill level of this soldier as a fraction of the maximum skill level
    float
    skill() const
    {
        SkillType max_skill = 0;
        max_skill = ~max_skill;
        return static_cast<float>(m_aggression) / max_skill;
    }

    unsigned char
    killRadius() const
    {
        return m_killRadiusMap[m_type];
    }

    // Destructor
    ~Soldier()
    { }

private:
    // creates soldier type to kill radius map
    static
    std::map<Type, unsigned char>
    createKillRadiusMap()
    {
        std::map<Type, unsigned char> killRadiusMap;
        killRadiusMap[LEADER] = 1;
        killRadiusMap[SWORDSMAN] = 1;
        killRadiusMap[ARCHER] = 5;
        return killRadiusMap;
    }

private:
    // stores kill radii for different soldier types
    static std::map<Type, unsigned char> m_killRadiusMap;

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

std::map<Soldier::Type, unsigned char> Soldier::m_killRadiusMap = Soldier::createKillRadiusMap();

std::ostream& operator<< (std::ostream& stream, const Soldier& t)
{
    if (t.empty())
        stream << " ";
    else if (t.army() == 0)
        stream << "x";
    else if (t.army() == 1)
        stream << "o";
    else
        stream << "-";
    return stream;
}

#endif // SOLDIER_H
