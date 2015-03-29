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

        // claim-a-cell loop
        for (std::size_t x = 0; x < m_nrows; ++x) {
            for (std::size_t y = 0; y < m_ncols; ++y) {
                // do the calculations only if there is a soldier in the current cell
                if (m_soldiers(x, y) != 0) {
                    calculateGlobalPreference(x, y, target_x, target_y, m_g);
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
                    // move the chosen soldier to this cell and add to the dynamic field in the previous cell

                    // now unset all the claimed bits for this cells
                    m_claimed(x, y) = 0;
                    // update the extended neighborhood counts
                }
            }
        }
    }

    void
    kill()
    {
        // choose-a-kill loop
        for (std::size_t x = 0; x < m_nrows; ++x) {
            for (std::size_t y = 0; y < m_ncols; ++y) {
            }
        }

        // actual kill loop
        for (std::size_t x = 0; x < m_nrows; ++x) {
            for (std::size_t y = 0; y < m_ncols; ++y) {
            }
        }
    }

    /// Destructor
    ~FloorField()
    { }

private:
    /// computes global preference matrix for a soldier, based on global target coordinates
    void
    calculateGlobalPreference(const std::size_t x, const std::size_t y,
                              const std::size_t target_x, const std::size_t target_y,
                              float (&m_g)[3][3]) const
    {
        for (unsigned char i = 0; i < 3; ++i) {
            for (unsigned char j = 0; j < 3; ++j) {
            }
        }
    }

    /// computes local preference matrix for a soldier, based on extended neighborhood
    void
    calculateLocalPreference(const std::size_t x, const std::size_t y, float (&m_l)[3][3]) const
    {
        for (unsigned char i = 0; i < 3; ++i) {
            for (unsigned char j = 0; j < 3; ++j) {
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
        // happiness of the soldier
        // TODO: calculate happiness based on extended neighborhood counts
        float h = 0.5;

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

        // store the probability with which the soldier wants to move to the target cell
        m_probability(x, y) = trans_prob[i][j];
        // also set the bit corresponding to this soldier in the target cell
        // first find the index of this cell relative to the target cell
        unsigned char index = 0;
        m_claimed(i, j) = m_claimed(i, j) | (1 << index);
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
