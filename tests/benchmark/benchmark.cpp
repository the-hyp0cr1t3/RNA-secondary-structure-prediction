#include <benchmark/benchmark.h>
#include <generators.hpp>

#undef RECOVER_MATCHING
#include <NucleicAcidSequence.hpp>
#include <Predictor.hpp>

static void BM_Random(benchmark::State& state) {
    int n = state.range(0);

    for(auto _ : state) {
        state.PauseTiming();

        RNA::NASeq seq = generators::gen_random(n);
        RNA::Predictor predictor(seq);

        state.ResumeTiming();

        size_t num_matchings = predictor.find_max_matching();
        benchmark::DoNotOptimize(num_matchings);
    }

    state.counters["n"] = n;
    state.SetComplexityN(n);
}

BENCHMARK(BM_Random)
    ->DenseRange(MIN_LEN, 500, 500 / 50)
    ->Complexity(benchmark::oNCubed);

BENCHMARK_MAIN();

/*
./bin/bench --benchmark_counters_tabular=true
./bin/bench --benchmark_counters_tabular=true --benchmark_format=csv
*/