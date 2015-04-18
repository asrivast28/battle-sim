/**
 * @file    floor_field.hpp
 * @author  Ankit Srivastava <asrivast@gatech.edu>
 * @brief
 *
 * Copyright (c) TODO
 */

#ifndef BATTLE_FIELD_H
#define BATTLE_FIELD_H

#include "matrix.hpp"
#include "soldier.hpp"

#include <cmath>
#include <limits>
#include <numeric>
#include <vector>


// DEBUG_MSG won't be printed in release build
#ifdef NDEBUG
#define DEBUG_MSG(format, ...)
#else
#define DEBUG_MSG(format, ...) fprintf(stderr, format, __VA_ARGS__)
#endif

/// returns a uniform random number in the range [0.0, 1.0]
/// TODO: replace with a random number generator which is actually uniform
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

// granulaty of measurement for float types
const float FLOAT_EPSILON = std::numeric_limits<float>::epsilon();
/// uniformly picks an index from an array consisting of probability distribution
size_t
pickIndex(float* const prob_dist, const size_t max_index)
{
    std::partial_sum(prob_dist, prob_dist + max_index, prob_dist);
    //if (std::fabs(cum_dist[max_index - 1] - 1.0) > FLOAT_EPSILON) {
        //throw std::runtime_error("Probability distribution doesn't add up to 1.0!");
    //}

    float prob = uniformRandom();
    float prev_val = 0.0;
    for (size_t i = 0; i < max_index; ++i) {
        if ((prob > prev_val) && ((prob <= prob_dist[i]))) {
            return i;
        }
        prev_val = prob_dist[i];
    }
    throw std::runtime_error("Something went wrong while picking index!");
}


class BattleField {
private:
  // extended neighborhood size
  static const unsigned char m_k = 5;
  // dynamic field decay constant
  static const unsigned char m_beta = 1;

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
    BattleField(const size_t nrows, const size_t ncols)
        : m_nrows(nrows), m_ncols(ncols), m_target(ANNIHILATE_ENEMY),
          m_status(ONGOING), m_winner(-1),
          m_soldiers(nrows, ncols), m_static(nrows, ncols, 255),
          m_claimed(nrows, ncols), m_probability(nrows, ncols),
          m_lastmove(nrows, ncols)
    {
        // initialize the target coordinates
        m_flag_x[0] = 0;
        m_flag_y[0] = 0;

        m_flag_x[1] = 0;
        m_flag_y[1] = 0;

        // initialize the dynamic field matrices
        m_dynamic[0] = matrix<unsigned char>(nrows, ncols, 1);
        m_dynamic[1] = matrix<unsigned char>(nrows, ncols, 1);

        // initialize the extended neighborhood counts
        m_neighbors[0] = matrix<unsigned char>(nrows, ncols);
        m_neighbors[1] = matrix<unsigned char>(nrows, ncols);

        m_total_soldiers[0] = 0;
        m_total_soldiers[1] = 0;
    }

    BattleField(size_t nrows, size_t ncols, unsigned char* accessibility)
        : m_nrows(nrows), m_ncols(ncols), m_target(ANNIHILATE_ENEMY),
          m_status(ONGOING), m_winner(-1),
          m_soldiers(nrows, ncols), m_static(nrows, ncols, accessibility, accessibility + nrows * ncols),
          m_claimed(nrows, ncols), m_probability(nrows, ncols), m_lastmove(nrows, ncols)
    {
        // initialize the target coordinates
        m_flag_x[0] = 0;
        m_flag_y[0] = 0;

        m_flag_x[1] = 0;
        m_flag_y[1] = 0;

        // initialize the dynamic field matrices
        m_dynamic[0] = matrix<unsigned char>(nrows, ncols, 1);
        m_dynamic[1] = matrix<unsigned char>(nrows, ncols, 1);

        // initialize the extended neighborhood counts
        m_neighbors[0] = matrix<unsigned char>(nrows, ncols);
        m_neighbors[1] = matrix<unsigned char>(nrows, ncols);

        m_total_soldiers[0] = 0;
        m_total_soldiers[1] = 0;
    }

    // access to the matrix of sodiers
    const matrix<Soldier>&
    mat() const
    {
        return m_soldiers;
    }

    matrix<Soldier>&
    mat()
    {
        return m_soldiers;
    }

    // print the matrix
    void
    printGrid() const {
        std::cout << m_soldiers << std::endl;
    }

    void
    setSoldiers(const std::vector<std::pair<size_t, Soldier> >& soldiers)
    {
        for (std::vector<std::pair<size_t, Soldier> >::const_iterator s = soldiers.begin(); s != soldiers.end(); ++s) {
            m_soldiers.at(s->first / m_ncols, s->first % m_ncols) = s->second;
        }
        initializeNeighborCounts();
    }

    void
    setTarget(const Target target)
    {
        m_target = target;
    }

    /// initializes the extended neighborhood counts
    void
    initializeNeighborCounts()
    {
        // calculate the initial count of soldiers of both the armies in k-neighborhood
        for (size_t x = 0; x < m_nrows; ++x) {
            for (size_t y = 0; y < m_ncols; ++y) {
                if (!m_soldiers(x, y).empty()) {
                    updateNeighborCounts(x, y, m_soldiers(x, y).army(), true);
                    ++m_total_soldiers[m_soldiers(x, y).army()];
                }
            }
        }
    }

    void
    setFlag(const unsigned char army, const size_t flag_x, const size_t flag_y)
    {
        m_flag_x[army] = flag_x;
        m_flag_y[army] = flag_y;
    }

    void
    move()
    {
        // global matrix of preference
        float mat_g[3 * 3] = {};
        // local matrix of preference
        float mat_l[3 * 3] = {};
        //last move matrix
        float mat_m[3 * 3] = {};
        // transitional probability matrix
        float trans_prob[3 * 3] = {};

        // claim-a-cell loop
        for (size_t x = 0; x < m_nrows; ++x) {
            for (size_t y = 0; y < m_ncols; ++y) {
                // do the calculations only if there is a soldier in the current cell
                if (!m_soldiers(x, y).empty()) {
                    // if the fourth param is true, goal is to destroy enemy. If nothing provided, false by default
                    calculateGlobalPreference(x, y, mat_g);
                    calculateLocalPreference(x, y, mat_l);
                    calculateMovementMatrix(x,y,mat_m);
                    calculateTransitionalProbabilities(x, y, mat_g, mat_l, mat_m, trans_prob);
                    claimCell(x, y, trans_prob);
                }
            }
        }

        // conflict resolution loop
        for (size_t x = 0; x < m_nrows; ++x) {
            for (size_t y = 0; y < m_ncols; ++y) {
                unsigned char claimed = m_claimed(x, y);
                // resolve conflict if this cell is claimed by more than one soldier
                if ((claimed > 0) && ((claimed & (claimed - 1)) != 0)) {
                    unsigned char max_size = sizeof(unsigned char) * 8;
                    float sum_prob = 0.0;
                    float rel_prob[max_size];
                    for (unsigned char a = 0, b = 1; a < max_size; ++a, b = b << 1) {
                        if ((claimed & b) != 0) {
                            // find the relative index of the cell that set this bit
                            unsigned char t = a + a / 4;
                            unsigned char i = 2 - (t / 3);
                            unsigned char j = 2 - (t % 3);
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
                // decay the dynamic field in this cell
                // minimum dynamic field in any cell is 1
                // we don't want dynamic field going to 0 now, do we?
                m_dynamic[0](x, y) = m_dynamic[0](x, y) * exp(-1.0 * static_cast<double>(m_beta));
                if (m_dynamic[0](x, y) == 0) {
                    m_dynamic[0](x, y) = 1;
                }
                m_dynamic[1](x, y) = m_dynamic[1](x, y) * exp(-1.0 * static_cast<double>(m_beta));
                if (m_dynamic[1](x, y) == 0) {
                    m_dynamic[1](x, y) = 1;
                }
            }
        }

        // actual movement loop
        for (size_t x = 0; x < m_nrows; ++x) {
            for (size_t y = 0; y < m_ncols; ++y) {
                // check if any soldier is moving to this cell
                if (m_claimed(x, y) > 0) {
                    // move the chosen soldier to this cell
                    unsigned char a = static_cast<unsigned char>(log2(m_claimed(x, y)));
                    unsigned char t = a + a / 4;
                    unsigned char i = 2 - (t / 3);
                    unsigned char j = 2 - (t % 3);

                    const Soldier& s = m_soldiers(x + i - 1, y + j - 1);
                    m_soldiers(x, y) = s;
                    m_soldiers(x + i - 1, y + j - 1).clear();
                    // update the move direction: 3*(t/3) + t%3 - (3*(t/3) + t%3)/5, translates to be 'a'
                    m_lastmove(x,y) = a;
		
                    // increase neighbor count in the new neighborhood
                    updateNeighborCounts(x, y, s.army(), true);
                    // decrease neighbor count in the old neighborhood
                    updateNeighborCounts(x + i - 1, y + j - 1, s.army(), false);

                    //DEBUG_MSG("Soldier at index (%zd, %zd) moving to (%zd, %zd)\n", x + i - 1, y + j - 1, x, y);

                    // add the dynamic field of this soldier in the previous cell
                    m_dynamic[s.army()](x + i - 1, y + j - 1) += s.field();

                    // now unset all the claimed bits for this cells
                    m_claimed(x, y) = 0;
                }
                // set probability to 0.0 for all the cells
                m_probability(x, y) = 0.0;
            }
        }
        // check if target has been achieved
        if ((m_target == CAPTURE_FLAG) && (m_status == ONGOING)) {
            bool captured[2] = {false, false};
            for (unsigned char i = 0; i < 2; ++i) {
                const Soldier& s = m_soldiers(m_flag_x[i], m_flag_y[i]);
                captured[i] = !s.empty() && (s.army() == i);
            }
            if (captured[0] && captured[1]) {
                m_status = TIED;
            }
            else if (captured[0] || captured[1]) {
                m_status = WON;
                m_winner = captured[0] ? 0 : 1;
            }
        }
    }

    void
    move(std::vector<std::pair<size_t, unsigned char> >& positions)
    {
        // first call actual move
        move();

        // now record positions
        size_t count = 0;
        for (size_t x = 0; x < m_nrows; ++x) {
            for (size_t y = 0; y < m_ncols; ++y) {
                if (!m_soldiers(x, y).empty()) {
                    const Soldier& s = m_soldiers(x, y);
                    unsigned char info = s.army() << 7;
                    info = info | static_cast<unsigned char>(s.type());
                    positions[count++] = std::make_pair(x * m_ncols + y, info);
                }
            }
        }
    }

    size_t
    kill(std::vector<size_t>& positions, bool record = true)
    {
        // choose-a-kill loop
        for (size_t x = 0; x < m_nrows; ++x) {
            for (size_t y = 0; y < m_ncols; ++y) {
                if (!m_soldiers(x, y).empty()) {
                    unsigned char army = m_soldiers(x, y).army();
                    unsigned char r = m_soldiers(x, y).killRadius();
                    std::vector<std::pair<unsigned char, unsigned char> > potentials;
                    // scan the cells in kill radius and record all the enemy soldiers
                    for (unsigned char i = 0; i < 2 * r; ++i) {
                        if ((x + i >= r) && (x + i - r < m_nrows)) {
                            for (unsigned char j = 0; j < 2 * r; ++j) {
                                if ((y + j >= r) && (y + j - r < m_ncols)) {
                                    if (m_soldiers(x + i - r, y + j - r).army() != army) {
                                        potentials.push_back(std::make_pair(i, j));
                                    }
                                }
                            }
                        }
                    }
                    if (potentials.size() > 0) {
                        // relative index of the chosen enemy
                        unsigned char index = static_cast<unsigned char>(randomAtMax(potentials.size() - 1));
                        unsigned char i = potentials[index].first;
                        unsigned char j = potentials[index].second;
                        // add to the kill probability of the soldier
                        m_probability(x + i - r, y + j - r) += m_soldiers(x, y).skill();
                    }
                }
            }
        }

        // actual kill loop
        size_t count = 0;
        for (size_t x = 0; x < m_nrows; ++x) {
            for (size_t y = 0; y < m_ncols; ++y) {
                if (!m_soldiers(x, y).empty() && (std::fabs(m_probability(x, y) - 0.0) > FLOAT_EPSILON)) {
                    float s_self = m_soldiers(x, y).skill();
                    float p_survival = s_self / (m_probability(x, y) + s_self);
                    // pick a uniform random number in the range [0.0, 1.0]
                    // kill the soldier if the chosen number is greater than the probability of survival
                    bool survives = uniformRandom() < p_survival;
                    if (!survives) {
                        // TODO: collect statistics for the soldier before killing
                        // kill the soldier
                        //DEBUG_MSG("killing (%zd, %zd)\n", x, y);
                        unsigned char army = m_soldiers(x, y).army();
                        m_soldiers(x, y).kill();
                        updateNeighborCounts(x, y, army, false);
                        if (record) {
                            if (count < positions.size()) {
                                positions[count++] = x * m_ncols + y;
                            }
                            else {
                                positions.push_back(x * m_ncols + y);
                                ++count;
                            }
                        }
                        --m_total_soldiers[army];
                    }
                    m_probability(x, y) = 0.0;
                }
            }
        }

        // check if target has been achieved
        if ((m_target == ANNIHILATE_ENEMY) && (m_status == ONGOING)) {
            if ((m_total_soldiers[0] == 0) && (m_total_soldiers[1] == 0)) {
                m_status = TIED;
            }
            else if ((m_total_soldiers[0] == 0) || (m_total_soldiers[1] == 0)) {
                m_status = WON;
                m_winner = (m_total_soldiers[0] == 0) ? 1 : 0;
            }
        }
        return count;
    }

    size_t
    kill()
    {
        std::vector<size_t> dummy;
        return kill(dummy, false);
    }

    size_t
    soldierCount(const unsigned char army) const
    {
        return m_total_soldiers[army];
    }

    Status
    status() const
    {
        return m_status;
    }

    unsigned char
    winner() const
    {
        return m_winner;
    }

    /// Destructor
    ~BattleField()
    { }

private:
    /// increase or decrease extended neighborhood counts for a cell
    void
    updateNeighborCounts(const size_t x, const size_t y, const unsigned char army, const bool increase)
    {
        for (unsigned char i = 0; i < 2 * m_k; ++i) {
            // decrease neighbor count in the neighborhood
            if ((x + i >= m_k) && (x + i - m_k < m_nrows)) {
                for (unsigned char j = 0; j < 2 * m_k; ++j) {
                    if ((y + j >= m_k) && (y + j - m_k < m_ncols)) {
                        // there is something wrong if we are trying to decrease neighbor count when it is 0
                        assert(increase || (m_neighbors[army](x + i - m_k, y + j - m_k) > 0));
                        m_neighbors[army](x + i - m_k, y + j - m_k) += (increase ? 1 : -1);
                    }
                }
            }
        }
    }

    /// computes global preference matrix for a soldier, based on global target coordinates
    void
    calculateGlobalPreference(const size_t x, const size_t y, float* const mat_g) const
    {
        // if we want to destroy the enemy, calculate in a different way
        if (m_target == ANNIHILATE_ENEMY) {
            // if there is no enemy nearby, false will be returned and the matrix will be calculated the usual way. Otherwise it will return here.
            bool should_use = calculateNearEnemyDirection(x,y,mat_g);
            if (should_use) {
                return;
            }
        }
	
        unsigned char army = m_soldiers(x, y).army();

        // give higher global preference to the cell which takes the soldier closer to the target
        float sum_distance = 0.0;
        float max_distance = 0.0;
        for (unsigned char i = 0; i < 3; ++i) {
              size_t diff_x = m_flag_x[army] > (x + i) ? (m_flag_x[army] - (x + i - 1)) : ((x + i - 1) - m_flag_x[army]);
              for (unsigned char j = 0; j < 3; ++j) {
                  size_t diff_y = m_flag_y[army] > (y + j) ? (m_flag_y[army] - (y + j - 1)) : ((y + j - 1) - m_flag_y[army]);
                  float distance = sqrt(pow(diff_x, 2.0) + pow(diff_y, 2.0));
                  mat_g[i * 3 + j] = distance;
                  sum_distance += distance;
                  if (max_distance < distance) {
                      max_distance = distance;
                  }
              }
        }

        // subtract the minimum distance from all the distances and then normalize
        sum_distance = (3 * 3 * max_distance) - sum_distance;
        if (std::fabs(sum_distance - 0.0) > FLOAT_EPSILON) {
            for (unsigned char i = 0; i < 3 * 3; ++i) {
                mat_g[i] = max_distance - mat_g[i];
                mat_g[i] /= sum_distance;
            }
        }
    }

    // support function for calculating global matrix of preference, in case we want to destroy opponent's army
    bool
    calculateNearEnemyDirection(const size_t x, const size_t y, float* const mat_g) const
    {
        // enemy army
        unsigned char e_army = m_soldiers(x,y).enemy();
        float sum_neighbors = 0.0;
        // first store the neighbor counts in the k-neighborhood
        for (unsigned char i = 0; i < 3; ++i) {
            // check the row bounds
            if ((x + i * m_k >= m_k) && (x +  i * m_k - m_k < m_nrows)) {
                for (unsigned char j = 0; j < 3; ++j) {
                    // check the column bounds
                    if ((y + j * m_k >= m_k) && (y + j * m_k - m_k < m_ncols)) {
                        mat_g[i * 3 + j] = m_neighbors[e_army](x + i * m_k - m_k, y + j * m_k - m_k);
                    }
                    sum_neighbors += mat_g[i * 3 + j];
                }
            }
        }

        // if there are no enemy neighbors nearby, should not be used
        if (sum_neighbors == 0) {
            return false;
        }
        if (std::fabs(sum_neighbors - 0.0) > FLOAT_EPSILON) {
            // now normalize based on the total sum of counts
            for (unsigned char i = 0; i < 3 * 3; ++i) {
                mat_g[i] /= sum_neighbors;
            }
        }
        return true;
    }

    /// computes local preference matrix for a soldier, based on extended neighborhood
    void
    calculateLocalPreference(const size_t x, const size_t y, float* const mat_l) const
    {
        unsigned char army = m_soldiers(x, y).army();

        float sum_neighbors = 0.0;
        // first store the neighbor counts in the k-neighborhood
        for (unsigned char i = 0; i < 3; ++i) {
            // check the row bounds
            if ((x + i * m_k >= m_k) && (x +  i * m_k - m_k < m_nrows)) {
                for (unsigned char j = 0; j < 3; ++j) {
                    // check the column bounds
                    if ((y + j * m_k >= m_k) && (y + j * m_k - m_k < m_ncols)) {
                        mat_l[i * 3 + j] = m_neighbors[army](x + i * m_k - m_k, y + j * m_k - m_k);
                    }
                    sum_neighbors += mat_l[i * 3 + j];
                }
            }
        }

        if (std::fabs(sum_neighbors - 0.0) > FLOAT_EPSILON) {
            // now normalize based on the total sum of counts
            for (unsigned char i = 0; i < 3 * 3; ++i) {
                mat_l[i] /= sum_neighbors;
            }
        }
    }

    void calculateMovementMatrix(const size_t x, const size_t y, float* const mat_m) const
    {
        unsigned char a = m_lastmove(x,y);
        unsigned char t = a + a/4;
        unsigned char i = t/3;
        unsigned char j = t%3;
        for (unsigned char p=0; p<3; p++) {
          for (unsigned char q=0; q<3; q++) {
            mat_m[p*3+q] = 0;
          }
        }
        mat_m[i*3+j]=1;
    }

    /// calculate transitional probabilities for a soldier, based on a soldier attributes, local and global preference matrix
    void
    calculateTransitionalProbabilities(const size_t x, const size_t y,
                                       float* const mat_g, float* const mat_l, float* const mat_m,
                                       float* const trans_prob) const
    {
        const Soldier& soldier = m_soldiers(x, y);
        // aggression of the soldier
        float a = soldier.aggression();
        // happiness of the soldier is defined as the relative count of soldiers of the same army
        float h = m_neighbors[soldier.army()](x, y);
        h /= (m_neighbors[soldier.army()](x, y) + m_neighbors[soldier.enemy()](x, y));
        float p=0.25;

        float sum_prob = 0.0;
        for (unsigned char i = 0; i < 3; ++i) {
            // check the row bounds
            if ((x + i >= 1) && (x + i - 1 < m_nrows)) {
                for (unsigned char j = 0; j < 3; ++j) {
                    // check the column bounds
                    if ((y + j >= 1) && (y + j - 1 < m_ncols)) {
                        if (m_soldiers(x + i - 1, y + j - 1).empty()) {
                            // calculate actual matrix of preference for this index
                            // TODO: refine the following expression to give preference to the current dir
                            float mat_ij = p * mat_m[i * 3 + j] + (1 - p) * (a * mat_g[i * 3 + j] + (1 - a) * (h * mat_g[i * 3 + j] + (1 - h) * mat_l[i * 3 + j]));
                            // calculate transitional probability
                            // TODO: refine the following expression?
                            assert (mat_ij >= 0);
                            trans_prob[i * 3 + j] = mat_ij * (m_dynamic[soldier.army()](x + i - 1, y + j - 1) / 255.0) * (m_static(x + i - 1, y + j - 1) / 255.0);

                            // reset global and local matrix of preference
                            mat_g[i * 3 + j] = 0.0;
                            mat_l[i * 3 + j] = 0.0;
                            mat_m[i * 3 + j] = 0.0;
                        }
                    }
                    sum_prob += trans_prob[i * 3 + j];
                }
            }
        }

        if (std::fabs(sum_prob - 0.0) < FLOAT_EPSILON) {
            // set probability of staying put to 1.0
            trans_prob[1 * 3 + 1] = 1.0;
        }
        else {
            // normalize the probabilities now so that the sum is 1.0
            for (unsigned char i = 0; i < 3 * 3; ++i) {
                trans_prob[i] /= sum_prob;
            }
        }
    }

    /// claim a cell for a soldier to which it wants to move in this timestep
    void
    claimCell(const size_t x, const size_t y, float* const trans_prob)
    {
        // relative index of the target cell
        unsigned char index = static_cast<unsigned char>(pickIndex(trans_prob, 3 * 3));

        // if the soldier wants to move to another cell, (i, j) != (1, 1), index = 4
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

        // reset transitional probabilities
        for (unsigned char i = 0; i < 3 * 3; ++i) {
            trans_prob[i] = 0.0;
        }
    }

private:
    // number of rows
    size_t m_nrows;
    // number of columns
    size_t m_ncols;

    // Type of target
    Target m_target;

    Status m_status;
    unsigned char m_winner;

    // flag coordinates
    size_t m_flag_x[2];
    size_t m_flag_y[2];

    // total number of alive soldiers
    size_t m_total_soldiers[2];

    // preference given to the last move of the soldier
    // float p;

    // grid for movement of the soldiers
    matrix<Soldier> m_soldiers;
    // matrix for storing static floor field
    matrix<unsigned char> m_static;
    // matrix for storing dynamic floor fields, for both the armies
    matrix<unsigned char> m_dynamic[2];

    // matrix for storing extended neighborhood counts, for both the armies
    matrix<unsigned char> m_neighbors[2];
    // matrix for storing soldiers which want to move to each cell
    matrix<unsigned char> m_claimed;
    // matrix for storing probabilities (kill and movement)
    matrix<float> m_probability;
    // matrix to store the last move
    matrix<unsigned char> m_lastmove;

}; // class BattleField

#endif // BATTLE_FIELD_H
