#pragma once

#include <NucleicAcidSequence.hpp>
#include <string>
#include <vector>

#define MIN_LEN 4

namespace RNA {

class Predictor {
    size_t n;
    NASeq naseq;

public:

    Predictor(const RNA::NASeq &seq);

    size_t find_max_matching();

    std::vector<std::pair<size_t, size_t>> recover_matching();

private:

    /// \f $dp_{i, j} \f $ = maximum number of base pairs in a secondary structure for \f $b_i b_{i+1} \ldots b_j \f $.
    std::vector<std::vector<size_t>> dp;

    std::vector<std::vector<int>> choices;

    void update_state(size_t l, size_t r, size_t m, size_t value);

};

} // RNA