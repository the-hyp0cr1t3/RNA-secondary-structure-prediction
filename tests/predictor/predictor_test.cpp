#include <gtest/gtest.h>
#include <NucleicAcidSequence.hpp>
#include <Predictor.hpp>
#include <algorithm>

namespace {

// ensures the matching does not break any rules ---------------------------------------------------
#define VALIDATE_MATCHING(matching, seq_len)                                                        \
    std::vector<int> frequency(seq_len);                                                            \
    int msz = matching.size();                                                                      \
    for(int i = 0; i < msz; i++) {                                                                  \
        auto [li, ri] = matching[i];                                                                \
                                                                                                    \
        /* bounds checking */                                                                       \
        EXPECT_LE(0, li);       /* 0 <= li */                                                       \
        EXPECT_LT(li, ri);      /* li < ri */                                                       \
        EXPECT_LT(ri, seq_len); /* ri < seq_len */                                                  \
                                                                                                    \
        /* no sharp turns: ends of each pair are separated by at least MIN_LEN intervening bases */ \
        EXPECT_LT(li, ri - MIN_LEN);                                                                \
                                                                                                    \
        /* elements in each pair consist of either {A,U} or {C,G} */                                \
        EXPECT_TRUE(RNA::matches(naseq[li], naseq[ri]));                                            \
                                                                                                    \
        frequency[li]++;                                                                            \
        frequency[ri]++;                                                                            \
                                                                                                    \
        for(int j = i + 1; j < msz; j++) {                                                          \
            auto [lj, rj] = matching[j];                                                            \
                                                                                                    \
            /* no knots: if (a,b) and (c,d) are two pairs, then we cannot have a < c < b < d */     \
            EXPECT_FALSE(li < lj and lj < ri and ri < rj                                            \
                            or lj < li and li < rj and rj < ri);                                    \
        }                                                                                           \
    }                                                                                               \
                                                                                                    \
    /* matching: no base appears in more than one pair */                                           \
    EXPECT_LE(*std::max_element(frequency.begin(), frequency.end()), 1);

// --------------------------------------------------------------------------------------------------


// macro to verify the output of a sequence ------------
#define DO_TEST(sequence, expected_num_matchings)       \
    RNA::NASeq naseq(sequence);                         \
    RNA::Predictor predictor(naseq);                    \
                                                        \
    int num_matchings = predictor.find_max_matching();  \
    auto matching = predictor.recover_matchings();      \
                                                        \
    EXPECT_EQ(num_matchings, expected_num_matchings);   \
    VALIDATE_MATCHING(matching, naseq.length())

// ------------------------------------------------------


TEST(PredictorTests, NDB7EFG) {
    DO_TEST("GGCGAAGAACCGGGGAGCC", 4);
}

TEST(PredictorTests, RandomLong) {
    DO_TEST("GGGUGAUUAGCUCAGCUGGGAGAGCACCUCCCUUACAAGGAGGGGGUCGGCGGUUCGAUCCCGUCAUCACCCACCA", 26);
}

TEST(PredictorTests, NDB6SVS) {
    DO_TEST("GAAGGUUUUUCUUUUUCCUGAAGGCGAAAGUCUUCAGGUUUUUCUUUUUGGCCUUUCUUAAAAAAAAAAAAAGAAAA", 28);
}

} // namespace