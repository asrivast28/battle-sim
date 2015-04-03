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

#include <numeric>
#include <vector>

/// returns a uniform random number in the range [0.0, 1.0]
/// TODO: replace with a more reliable random number generator
float
uniformRandom()
{
    return rand() / static_cast<float>(RAND_MAX);
}

// slightly modified version of http://stackoverflow.com/a/6852396
// Assumes 0 <= max <= RAND_MAX
// Returns in the half-open interval [0, max]
int
randomAtMax(int max) {
    unsigned num_bins = static_cast<unsigned>(max) + 1;
    // max <= RAND_MAX < ULONG_MAX, so this is okay.
    unsigned num_rand = static_cast<unsigned>(RAND_MAX) + 1;
    unsigned bin_size = num_rand / num_bins;
    unsigned defect = num_rand % num_bins;

    int x;
    do {
        x = rand();
    }
    // This is carefully written not to overflow
    while ((num_rand - defect) <= static_cast<unsigned>(x));

    // Truncated division is intentional
    return x / static_cast<int>(bin_size);
}



/// uniformly picks an index from an array consisting of probability distribution
std::size_t
pickIndex(float* const prob_dist, const std::size_t max_index)
{
    std::partial_sum(prob_dist, prob_dist + max_index, prob_dist);
    if (prob_dist[max_index - 1] != 1.0) {
        throw std::runtime_error("Probability distribution doesn't add up to 1.0!");
    }
    double prob = uniformRandom();
    double prev_val = 0.0;
    for (std::size_t i = 0; i < max_index; ++i) {
        if ((prob > prev_val) && (prob <= prob_dist[i])) {
            return i;
        }
    }
    throw std::runtime_error("Something went wrong while picking index!");
}


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

    // access to the matrix of sodiers
    const matrix<Soldier>&
    mat() const {
      return m_soldiers;
    }

    matrix<Soldier>&
    mat() {
      return m_soldiers;
    }

    void
    move()
    {
        // global matrix of preference
        float m_g[3 * 3];
        // local matrix of preference
        float m_l[3 * 3];
        // transitional probability matrix
        float trans_prob[3 * 3];

        // calculate the count of soldiers of both the armies in k-neighborhood
        // TODO: do this only once in the beginning and then only update later
        for (std::size_t x = 0; x < m_nrows; ++x) {
            for (std::size_t y = 0; y < m_ncols; ++y) {
                unsigned char neighbors[] = {0, 0};
                for (unsigned char i = 0; i < 2 * m_k; ++i) {
                    if ((x + i >= m_k) && (x + i - m_k < m_nrows)) {
                        for (unsigned char j = 0; j < 2 * m_k; ++j) {
                            if ((y + j >= m_k) && (y + j - m_k < m_ncols)) {
                                ++neighbors[m_soldiers(x + i - m_k, y + j - m_k).army()];
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
                if (!m_soldiers(x, y).empty()) {
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
                unsigned char claimed = m_claimed(x, y);
                // resolve conflict only if this cell is claimed by more than one soldiers
                if ((claimed > 0) && ((claimed & (claimed - 1)) != 0)) {
                    unsigned char max_size = sizeof(unsigned char) * 8;
                    float sum_prob = 0.0;
                    float rel_prob[max_size];
                    for (unsigned char a = 0, b = 1; a < max_size; ++a, b = b << 1) {
                        if ((claimed & b) != 0) {
                            // find the relative index of this cell relative to the cell that set this bit
                            unsigned char s = a + a / 4;
                            unsigned char i = s / 3;
                            unsigned char j = s % 3;
                            // now calculate the relative index of the cell that set this bit
                            i = (i + 1) % 3;
                            j = (j + 1) % 3;
                            // store the probability with which this bit was set
                            rel_prob[a] = m_probability(x + i - 1, y + j - 1);
                        }
                        else {
                            rel_prob[a] = 0.0;
                        }
                        sum_prob += rel_prob[a];
                    }
                    // normalize the probabilities
                    for (unsigned char a = 0; a < 8; ++a) {
                        rel_prob[a] /= sum_prob;
                    }
                    // choose one of the claimants as the "lucky" one, based on the relative probabilities
                    unsigned char index = static_cast<unsigned char>(pickIndex(rel_prob, max_size));
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
                if (!m_soldiers(x, y).empty()) {
                    unsigned char army = m_soldiers(x, y).army();
                    unsigned char k = m_soldiers(x, y).killRadius();
                    std::vector<std::pair<unsigned char, unsigned char> > potentials;
                    // scan the cells on kill rectangle and record all the enemy soldiers
                    for (unsigned char i = 0; i < 2; ++i) {
                        if (x + (i * k) >= k) {
                            for (unsigned char j = 0; j < 2 * k; ++j) {
                                if (y + j >= k && (y + j - k < m_nrows)) {
                                    const Soldier& soldier = m_soldiers(x + (i * k) - k, y + j - k);
                                    if (soldier.army() != army) {
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
                                    const Soldier& soldier = m_soldiers(x + i - k, y + (j * k) - k);
                                    if (soldier.army() != army) {
                                        potentials.push_back(std::make_pair(i, j * k));
                                    }
                                }
                            }
                        }
                    }
                    // relative index of the chosen enemy
                    unsigned char index = static_cast<unsigned char>(randomAtMax(potentials.size()));
                    unsigned char i = potentials[index].first;
                    unsigned char j = potentials[index].second;
                    // add to the kill probability of the soldier
                    m_probability(x + i - k, y + j - k) += m_soldiers(x, y).skill();
                }
            }
        }

        // actual kill loop
        for (std::size_t x = 0; x < m_nrows; ++x) {
            for (std::size_t y = 0; y < m_ncols; ++y) {
                if (!m_soldiers(x, y).empty()) {
                    float s_self = m_soldiers(x, y).skill();
                    float p_survival = s_self / (m_probability(x, y) + s_self);
                    // pick a uniform random number in the range [0.0, 1.0]
                    // kill the soldier if the chosen number is greater than the probability of survival
                    bool survives = uniformRandom() < p_survival;
                    if (!survives) {
                        // TODO: collect statistics for the soldier before killing
                        // kill the soldier
                        m_soldiers(x, y).kill();
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
    calculateGlobalPreference(const std::size_t x, const std::size_t y, float* const m_g) const
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
    calculateLocalPreference(const std::size_t x, const std::size_t y, float* const m_l) const
    {
        unsigned char army = m_soldiers(x, y).army();

        float sum = 0.0;
        // first store the neighbor counts in the k-neighborhood
        for (unsigned char i = 0; i < 3; ++i) {
            for (unsigned char j = 0; j < 3; ++j) {
                if (((x + i) > m_k) && ((y + j) > m_k)) {
                    m_l[i * 3 + j] = m_neighbors[army](x + i - 1 + m_k, y + j - 1 + m_k);
                }
                else {
                    m_l[i * 3 + j] = 0.0;
                }
                sum += m_l[i * 3 + j];
            }
        }

        // now normalize based on the total sum of counts
        // TODO: ensure that all the values in the matrix are non-zero
        for (unsigned char i = 0; i < 3; ++i) {
            for (unsigned char j = 0; j < 3; ++j) {
                m_l[i * 3 + j] /= sum;
            }
        }
    }

    /// calculate transitional probabilities for a soldier, based on a soldier attributes, local and global preference matrix
    void
    calculateTransitionalProbabilities(const std::size_t x, const std::size_t y,
                                       const float* const m_g, const float* const m_l,
                                       float* const trans_prob) const
    {
        const Soldier& soldier = m_soldiers(x, y);
        // aggression of the soldier
        float a = soldier.aggression();
        // happiness of the soldier is defined as the relative count of soldiers of the same army
        float h = m_neighbors[soldier.army()](x, y);
        h /= (m_neighbors[soldier.army()](x, y) + m_neighbors[soldier.enemy()](x, y));

        float sum_prob = 0.0;
        for (unsigned char i = 0; i < 3; ++i) {
            for (unsigned char j = 0; j < 3; ++j) {
                // check if the target cell is within bounds and is empty
                if (((x + i) > 0) && ((y + j) > 0) && ((x + i - 1) < m_nrows) && ((y + j - 1) < m_ncols) && (!m_soldiers(x + i - 1,  y + j - 1).empty())) {
                    // calculate actual matrix of preference for this index
                    float m_ij = a * m_g[i * 3 + j] + (1 - a) * (h * m_g[i * 3 + j] + (1 - h) * m_l[i * 3 + j]);
                    // calculate transitional probability
                    // TODO: use an accurate expression for calculating the probability
                    trans_prob[i * 3 + j] = m_ij * m_dynamic[soldier.army()](x + i - 1, y + j - 1) * m_static(x + i - 1, y + j - 1);
                }
                else {
                    // set transitional probability to 0.0 if the target cell is out of bounds or is occupied
                    trans_prob[i * 3 + j] = 0.0;
                }
                sum_prob += trans_prob[i * 3 + j];
            }
        }

        // normalize the probabilities now so that the sum is 1.0
        for (unsigned char i = 0; i < 3 * 3; ++i) {
            trans_prob[i] /= sum_prob;
        }
    }

    /// claim a cell for a soldier to which it wants to move in this timestep
    void
    claimCell(const std::size_t x, const std::size_t y, float* const trans_prob)
    {
        // relative index of the target cell
        unsigned char index = static_cast<unsigned char>(pickIndex(trans_prob, 3 * 3));

        // if the soldier wants to move to another cell, (i, j) != (1, 1)
        if (index != 4) {
            unsigned char i = index / 3;
            unsigned char j = index % 3;
            // store the probability with which the soldier wants to move to the target cell
            m_probability(x, y) = trans_prob[i * 3 + j];
            // also set the bit corresponding to this soldier in the target cell
            // first find the index of this cell relative to the target cell
            unsigned char bit_index = (3 * i + j) - ((3 * i + j) / 5);
            m_claimed(x + i - 1, y + j - 1) = m_claimed(x + i - 1, y + j - 1) | (1 << bit_index);
        }
    }

private:
    // number of rows
    std::size_t m_nrows;
    // number of columns
    std::size_t m_ncols;

    // grid for movement of the soldiers
    matrix<Soldier> m_soldiers;
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
