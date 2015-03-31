/**
 * @file    floor_field.hpp
 * @author  Ankit Srivastava <asrivast@gatech.edu>
 * @brief
 *
 * Copyright (c) TODO
 */

#ifndef FLOOR_FIELD_H
#define FLOOR_FIELD_H

#include "matrix.hpp"
#include "soldier.hpp"

#include <vector>

class FloorField {
private:
  // extended neighborhood size
  static const unsigned char m_k = 5;

public:
    FloorField(std::size_t nrows, std::size_t ncols)
        : m_nrows(nrows), m_ncols(ncols),
          m_soldiers(nrows, ncols), m_static(nrows, ncols),
          m_claimed(nrows, ncols), m_probability(nrows, ncols)
    {
        // initialize the dynamic field matrices
        m_dynamic[0] = matrix<float>(nrows, ncols);
        m_dynamic[1] = matrix<float>(nrows, ncols);

        // initialize the extended neighborhood counts
        m_neighbors[0] = matrix<unsigned char>(nrows, ncols);
        m_neighbors[1] = matrix<unsigned char>(nrows, ncols);

    }

    void
    move()
    {
        // global matrix of preference
        float m_g[3][3];
        // local matrix of preference
        float m_l[3][3];
        // transitional probability matrix
        float trans_prob[3][3];

        // calculate the count of soldiers of both the armies in k-neighborhood
        // TODO: do this only once in the beginning and then only update later
        for (std::size_t x = 0; x < m_nrows; ++x) {
            for (std::size_t y = 0; y < m_ncols; ++y) {
                unsigned char neighbors[] = {0, 0};
                for (unsigned char i = 0; i < 2 * m_k; ++i) {
                    if (x + i >= m_k) {
                        for (unsigned char j = 0; j < 2 * m_k; ++j) {
                            if (y + j >= m_k) {
                                ++neighbors[m_soldiers(x + i - m_k, y + j - m_k)->army()];
                            }
                        }
                    }
                }
                m_neighbors[0](x, y) = neighbors[0];
                m_neighbors[1](x, y) = neighbors[1];
            }
        }

        // claim-a-cell loop
        for (std::size_t x = 0; x < m_nrows; ++x) {
            for (std::size_t y = 0; y < m_ncols; ++y) {
                // do the calculations only if there is a soldier in the current cell
                if (m_soldiers(x, y) != 0) {
                    calculateGlobalPreference(x, y, m_g);
                    calculateLocalPreference(x, y, m_l);
                    calculateTransitionalProbabilities(x, y, m_g, m_l, trans_prob);
                    claimCell(x, y, trans_prob);
                }
            }
        }

        // conflict resolution loop
        for (std::size_t x = 0; x < m_nrows; ++x) {
            for (std::size_t y = 0; y < m_ncols; ++y) {
                // check if this cell is claimed by more than one soldiers
                unsigned char claimed = m_claimed(x, y);
                if ((claimed > 0) && ((claimed & (claimed - 1)) != 0)) {
                    unsigned char rel_prob[sizeof(unsigned char) * 8];
                    unsigned char b = 1;
                    for (unsigned char a = 1; a <= 8; ++a) {
                        if ((claimed & b) != 0) {
                            unsigned char i = a / 3;
                            unsigned char j = a % 3;
                            rel_prob[a] = m_probability(x + i - 1, y + j - 1);
                        }
                        else {
                            rel_prob[a] = 0.0;
                        }
                        b = b << 1;
                    }
                    unsigned char index = 8;
                    // TODO: choose one of the claimants as the "lucky" one, based on the relative probabilities
                    m_claimed(x, y) = 1 << index;
                }
                // TODO: decay the dynamic field in this cell
            }
        }

        // actual movement loop
        for (std::size_t x = 0; x < m_nrows; ++x) {
            for (std::size_t y = 0; y < m_ncols; ++y) {
                // check if any soldier is moving to this cell
                if (m_claimed(x, y) > 0) {
                    // move the chosen soldier to this cell
                    // add to the dynamic field in the previous cell

                    // now unset all the claimed bits for this cells
                    m_claimed(x, y) = 0;
                    // update the extended neighborhood counts
                }
                // set probability to 0.0 for all the cells
                m_probability(x, y) = 0.0;
            }
        }
    }

    void
    kill()
    {
        // choose-a-kill loop
        for (std::size_t x = 0; x < m_nrows; ++x) {
            for (std::size_t y = 0; y < m_ncols; ++y) {
                if (m_soldiers(x, y) != 0) {
                    unsigned char army = m_soldiers(x, y)->army();
                    unsigned char k = m_soldiers(x, y)->killRadius();
                    std::vector<std::pair<unsigned char, unsigned char> > potentials;
                    // scan the cells on kill rectangle and record all the enemy soldiers
                    for (unsigned char i = 0; i < 2; ++i) {
                        if (x + (i * k) >= k) {
                            for (unsigned char j = 0; j < 2 * k; ++j) {
                                if (y + j >= k) {
                                    const Soldier* const soldier = m_soldiers(x + (i * k) - k, y + j - k);
                                    if (soldier->army() != army) {
                                        potentials.push_back(std::make_pair(i * k, j));
                                    }
                                }
                            }
                        }
                    }
                    for (unsigned char j = 0; j < 2; ++j) {
                        if (y + (j * k) >= k) {
                            for (unsigned char i = 0; i < 2 * k; ++i) {
                                if (x + i >= k) {
                                    const Soldier* const soldier = m_soldiers(x + i - k, y + (j * k) - k);
                                    if (soldier->army() != army) {
                                        potentials.push_back(std::make_pair(i, j * k));
                                    }
                                }
                            }
                        }
                    }
                    // relative index of the chosen enemy
                    unsigned char i = 0, j = 0;
                    // TODO: uniformly pick one soldier to attack out of all the enemies in kill range
                    // add to the kill probability of the soldier
                    m_probability(x + i - k, y + j - k) += m_soldiers(x, y)->skill();
                }
            }
        }

        // actual kill loop
        for (std::size_t x = 0; x < m_nrows; ++x) {
            for (std::size_t y = 0; y < m_ncols; ++y) {
                if (m_soldiers(x, y) != 0) {
                    float s_self = m_soldiers(x, y)->skill();
                    float p_survival = s_self / (m_probability(x, y) + s_self);
                    bool survives = true;
                    // TODO: decide if the soldier survives or not, based on the survival probability
                    if (!survives) {
                        // TODO: collect statistics for the soldier
                        // free the cell
                        delete m_soldiers(x, y);
                        m_soldiers(x, y) = 0;
                    }
                }
            }
        }
    }

    /// Destructor
    ~FloorField()
    { }

private:
    /// computes global preference matrix for a soldier, based on global target coordinates
    void
    calculateGlobalPreference(const std::size_t x, const std::size_t y, float (&m_g)[3][3]) const
    {
        // target coordinates
        std::size_t target_x = 0, target_y = 0;
        // TODO: calculate target coordinates
        for (unsigned char i = 0; i < 3; ++i) {
            for (unsigned char j = 0; j < 3; ++j) {
            }
        }
    }

    /// computes local preference matrix for a soldier, based on extended neighborhood
    void
    calculateLocalPreference(const std::size_t x, const std::size_t y, float (&m_l)[3][3]) const
    {
        unsigned char army = m_soldiers(x, y)->army();

        float sum = 0.0;
        // first store the neighbor counts in the k-neighborhood
        for (unsigned char i = 0; i < 3; ++i) {
            for (unsigned char j = 0; j < 3; ++j) {
                if (((x + i) > m_k) && ((y + j) > m_k)) {
                    m_l[i][j] = m_neighbors[army](x + i - 1 + m_k, y + j - 1 + m_k);
                }
                else {
                    m_l[i][j] = 0.0;
                }
                sum += m_l[i][j];
            }
        }

        // now normalize based on the total sum of counts
        // TODO: ensure that all the values in the matrix are non-zero
        for (unsigned char i = 0; i < 3; ++i) {
            for (unsigned char j = 0; j < 3; ++j) {
                m_l[i][j] /= sum;
            }
        }
    }

    /// calculate transitional probabilities for a soldier, based on a soldier attributes, local and global preference matrix
    void
    calculateTransitionalProbabilities(const std::size_t x, const std::size_t y,
                                       const float (&m_g)[3][3], const float (&m_l)[3][3],
                                       float (&trans_prob)[3][3]) const
    {
        const Soldier* const soldier = m_soldiers(x, y);
        // aggression of the soldier
        float a = soldier->aggression();
        // happiness of the soldier is defined as the relative count of soldiers of the same army
        float h = m_neighbors[soldier->army()](x, y);
        h /= (m_neighbors[soldier->army()](x, y) + m_neighbors[soldier->enemy()](x, y));

        float sum_prob = 0.0;
        for (unsigned char i = 0; i < 3; ++i) {
            for (unsigned char j = 0; j < 3; ++j) {
                // check if the target cell is within bounds and is empty
                if (((x + i) > 0) && ((y + j) > 0) && (m_soldiers(x + i - 1,  y + j - 1) == 0)) {
                    // calculate actual matrix of preference for this index
                    float m_ij = a * m_g[i][j] + (1 - a) * (h * m_g[i][j] + (1 - h) * m_l[i][j]);
                    // calculate transitional probability
                    // TODO: use an accurate expression for calculating the probability
                    trans_prob[i][j] = m_ij * m_dynamic[soldier->army()](x + i - 1, y + j - 1) * m_static(x + i - 1, y + j - 1);
                }
                else {
                    // set transitional probability to 0.0 if the target cell is out of bounds or is occupied
                    trans_prob[i][j] = 0.0;
                }
                sum_prob += trans_prob[i][j];
            }
        }

        // normalize the probabilities now so that the sum is 1.0
        for (unsigned char i = 0; i < 3; ++i) {
            for (unsigned char j = 0; j < 3; ++j) {
                trans_prob[i][j] /= sum_prob;
            }
        }
    }

    /// claim a cell for a soldier to which it wants to move in this timestep
    void
    claimCell(const std::size_t x, const std::size_t y, const float (&trans_prob)[3][3])
    {
        // relative index of the target cell
        unsigned char i = 0, j = 0;
        // TODO: pick a cell (i, j) to move to, based on transitional probabilities

        // if the soldier wants to move to another cell, (i, j) != (1, 1)
        if ((i != 1) || (j != 1)) {
            // store the probability with which the soldier wants to move to the target cell
            m_probability(x, y) = trans_prob[i][j];
            // also set the bit corresponding to this soldier in the target cell
            // first find the index of this cell relative to the target cell
            unsigned char index = (3 * i + j) - ((3 * i + j) / 5);
            m_claimed(x + i - 1, y + j - 1) = m_claimed(x + i - 1, y + j - 1) | (1 << index);
        }
    }

private:
    // number of rows
    std::size_t m_nrows;
    // number of columns
    std::size_t m_ncols;

    // grid for movement of the soldiers
    matrix<Soldier*> m_soldiers;
    // matrix for storing static floor field
    matrix<float> m_static;
    // matrix for storing dynamic floor fields, for both the armies
    matrix<float> m_dynamic[2];

    // matrix for storing extended neighborhood counts, for both the armies
    matrix<unsigned char> m_neighbors[2];
    // matrix for storing soldiers which want to move to each cell
    matrix<unsigned char> m_claimed;
    // matrix for storing probabilities (kill and movement)
    matrix<float> m_probability;

}; // class FloorField

#endif // FLOOR_FIELD_H
