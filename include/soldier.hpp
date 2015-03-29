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
    /// default constructor
    Soldier()
      : m_army(), m_type(),
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
