#include <Predictor.hpp>
#include <stack>
#include <cassert>

RNA::Predictor::Predictor(const RNA::NASeq &seq)
    : n(seq.length()), naseq(seq),
        dp(seq.length(), std::vector<size_t>(seq.length())),
        choices(seq.length(), std::vector<int>(seq.length(), -1)) {}

size_t RNA::Predictor::find_max_matching() {

    for(size_t len = MIN_LEN + 1; len < n; len++) {
        for(size_t l = 0, r = len; r < n; l++, r++) {
            update_state(l, r, r - 1, dp[l][r - 1]);
            for(size_t m = l; m + MIN_LEN < r; m++) {
                size_t nval = (m > 0? dp[l][m - 1] : 0) + matches(naseq[m], naseq[r]) + dp[m + 1][r - 1];
                update_state(l, r, m, nval);
            }
        }
    }

    return dp[0][n - 1];
}

void RNA::Predictor::update_state(size_t l, size_t r, size_t m, size_t value) {

    if(value > dp[l][r]) {
        dp[l][r] = value;

#ifdef RECOVER_MATCHING
        choices[l][r] = m;
#endif

    }

}

std::vector<std::pair<size_t, size_t>> RNA::Predictor::recover_matchings() {
#ifndef RECOVER_MATCHING
    throw std::runtime_error("Attempt to recover choices when RECOVER_MATCHING is not defined");
#endif

    std::vector<std::pair<size_t, size_t>> matchings;
    matchings.reserve(dp[0][n - 1]);

    std::stack<std::pair<int, int>> pending {{ { 0, n - 1 } }};

    while(!pending.empty()) {
        auto [l, r] = pending.top();
        pending.pop();

        if(l < 0 or r < l or (int)n <= r) continue;

        int m = choices[l][r];

        if(m == -1) continue;

        if(m == (int)r - 1) {
            pending.emplace(l, r - 1);

        } else {
            pending.emplace(l, m - 1);
            pending.emplace(m + 1, r - 1);
            matchings.emplace_back(m, r);
        }
    }

    assert(matchings.size() == dp[0][n - 1]);

    return matchings;
}
