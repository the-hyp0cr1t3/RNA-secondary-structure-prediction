/**
 * @file Predictor.hpp
 * @author the-hyp0cr1t3
 * @brief Describes the predictor for secondary structure of RNA
 * @date 2022-04-18
 */
#pragma once

#include <NucleicAcidSequence.hpp>
#include <string>
#include <vector>

/// The minimum number of intervening bases between a pair in the matching
#define MIN_LEN 4

namespace RNA {

class Predictor {
    size_t n;       ///< The length of the nucleic acid base sequence
    NASeq naseq;    ///< The nucleic acid base sequence

public:

    /**
     * @brief Constructor
     *
     * @param seq The nucleic acid base sequence
     */
    Predictor(const RNA::NASeq &seq);

    /**
     * @brief Finds the cardinality of the maximum matching of base pairs
     *
     * Employs (iterative) dynamic programming on intervals to find
     * the cardinality of the maximum matching of base pairs that satisfy the following rules
     * * **No sharp turns:** The ends of each pair are separated by at least `MIN_LEN` intervening bases
     *   i.e. if \f$(i, j) \in S \f$, then \f$ i \lt j - \texttt{MIN_LEN} \f$.
     * * **Complementary base pairs:** The elements in each pair in S consist of either \f$ \{A,U\} \f$ or \f$ \{C,G\} \f$ (in either order).
     * * **Matching:** No base appears in more than one pair.
     * * **No knots:** If \f$(i,j) \in S\f$ and \f$(k,l) \in S\f$, then we cannot have \f$i \lt k \lt j \lt l\f$.
     *
     * \f$ dp_{\, l, \, r } \f$ denotes the maximum number of disjoint base pairs in the interval \f$ seq_{\, l \ldots r} \f$.
     *
     * \f$
     * dp_{\, l, \, r } = \max
     * \begin{cases}
     *     \, \, \, \, \, \, \, \, \, \, \, \, \, \, \, \, \, \, \, \, \, \, \,
     *     \, \, \, \, \, \, \, \, \, \, \, \, \, \, 0  &  \, r - l \le \texttt{MIN_LEN} \\
     *     \, \, \, \, \, \, \, \, \, \, \, \, \, \, \, \, \, \, \, \,
     *     \, \, \, \, \, \, \, \, \, \, \, \, \, \, dp_{\, l, \, r-1 }  &  \, r - l \gt \texttt{MIN_LEN} \\
     *     \underset{l \le m \le r}{\max} \, \, \, dp_{\, l, \, m-1 } + 1 + dp_{\, m+1, \, r-1}
     *      &  \, b_m \text{ and } b_r \text{ are complementary bases}
     * \end{cases}
     * \f$
     *
     * The time complexity of the algorithm is \f$ \mathcal{O}(n^3)\f$
     *
     * @return `size_t` The number of matching base pairs, i.e. \f$ dp_{\, l, \, r } \f$.
     */
    size_t find_max_matching();

    /**
     * @brief Recovers the indices of the matching base pairs
     * @pre `RECOVER_MATCHING` must be defined
     * @pre Must be called after `find_max_matching()`
     *
     * Every time a dp state is updated in `find_max_matching()`,
     * the choices table is updated with the state that made the transition. <br>
     * By backtracking through choices (previous states) starting from the final state \f$ dp_{\, 1, \, n } \f$,
     * all base pairs from the maximum matching can be recovered.
     *
     * \f$ choices_{\, l, \, r } \f$ denotes the choice (previous state) that led to \f$  dp_{\, l, \, r } \f$.
     *
     * \f$
     * choices_{\, l, \, r } = \max
     * \begin{cases}
     *     \, \, \, m &  \, dp_{\, l, \, m-1 } + 1 + dp_{\, m+1, \, r-1} \text{   is max} \\
     *     r-1  &  \, dp_{\, l, \, r-1 } \text{   is max} \\
     *     \, \, -1 &  \, \text{otherwise}
     * \end{cases}
     * \f$
     *
     * @return `std::vector<std::pair<size_t, size_t>>` The indices of the matching base pairs.
     */
    std::vector<std::pair<size_t, size_t>> recover_matchings();

private:

    /**
     * @brief The dp table which keeps track of optimal values for each state
     *
     * \f$ dp_{\, l, \, r } \f$ denotes the maximum number of disjoint base pairs in the interval \f$ seq_{\, l \ldots r} \f$.
     */
    std::vector<std::vector<size_t>> dp;

    /**
     * @brief The choices table which keeps track of the choices made (previous state) after every successful transition
     *
     * \f$ choices_{\, l, \, r } \f$ denotes the choice (previous state) that led to \f$  dp_{\, l, \, r } \f$.
     */
    std::vector<std::vector<int>> choices;

    /**
     * @brief Updates the dp table and conditionally the choices table if `RECOVER_MATCHING`  is defined
     *
     * @param l The left endpoint of the interval
     * @param r The right endpoint of the interval
     * @param m An interior point in the interval corresponding to \f$ choices_{\, l, \, r } \f$
     * @param value A candidate for the new value of \f$ dp_{\, l, \, r } \f$
     */
    void update_state(size_t l, size_t r, size_t m, size_t value);

};

} // RNA