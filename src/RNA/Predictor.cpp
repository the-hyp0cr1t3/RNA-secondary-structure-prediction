#include <Predictor.hpp>
#include <cassert>

RNA::Predictor::Predictor(const RNA::NASeq &seq): n(seq.length()), naseq(seq) {}

size_t RNA::Predictor::find_max_matching() {
    for(size_t i = 0; i < n; i++)
        for(size_t j = 0; j < n; j++)
            dp[i][j] = 0, choices[i][j] = -1;

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

RNA::pair *RNA::Predictor::recover_matchings() {
#ifndef RECOVER_MATCHING
    throw std::runtime_error("Attempt to recover choices when RECOVER_MATCHING is not defined");
#endif

    size_t idx = 0;
    RNA::pair *matchings = new RNA::pair[dp[0][n - 1]];

    int p_top = -1;
    RNA::pair pending[dp[0][n - 1]];

    auto push = [&pending, &p_top](int x, int y) {
        pending[++p_top] = { x, y };
    };

    push(0, n - 1);

    while(p_top != -1) {
        int l = pending[p_top].x, r = pending[p_top].y;
        --p_top;

        if(l < 0 or r < l or (int)n <= r) continue;

        int m = choices[l][r];

        if(m == -1) continue;

        if(m == (int)r - 1) {
            push(l, r - 1);

        } else {
            push(l, m - 1);
            push(m + 1, r - 1);
            matchings[idx++] = { m, r };
        }
    }

    assert(idx == dp[0][n - 1]);

    return matchings;
}
