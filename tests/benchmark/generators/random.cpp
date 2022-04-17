#include <generators.hpp>
#include <chrono>
#include <random>
#include <algorithm>

// random engine
std::mt19937 rng(std::chrono::steady_clock::now().time_since_epoch().count());

// [L, R] inclusive
auto rand_int = [](int L, int R) {
    return std::uniform_int_distribution<int>(L, R)(rng);
};

RNA::NASeq generators::gen_random(size_t len) {
    std::string s(len, 0);
    std::generate(s.begin(), s.end(), [] { return "AUCG"[rand_int(0, 3)]; });

    return RNA::NASeq(s);
}
