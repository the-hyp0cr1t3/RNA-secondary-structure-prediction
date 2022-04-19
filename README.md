# RNA secondary structure prediction
**Problem:** Given a nucleic acid sequence of RNA, find a [maximum matching](https://en.wikipedia.org/wiki/Matching_(graph_theory)) of $\{A,U\}$ or $\{C,G\} $ base pairs without knots or sharp turns.

This is a modern C++ implementation that employs (iterative) [dynamic programming](https://en.wikipedia.org/wiki/Dynamic_programming) on intervals to find the cardinality of the maximum matching of base pairs as well as the base pairs in the matching.

View the report [here](https://lucent-lebkuchen-97223f.netlify.app/report.html).

## Install Dependencies
To build the project you must have [CMake](https://cmake.org/install/) installed.

To install python dependencies
```sh
pip install -r requirements.txt
```

## Build Project
#### Configure:
```cmake
cmake -DCMAKE_BUILD_TYPE=Release -S . -B build
```
#### Build:
```cmake
cmake --build build --config Release
```

## Build Documentation
To build the documentation you must have [Doxygen](https://github.com/doxygen/doxygen) installed.

```cmake
cmake --build build --target docs
```

Output will be in `docs/html`.

## Usage

```bash
./bin/app [inputfile]
```

#### Example:
```
./bin/app sample.txt
```

**Note:** `inputfile` may also be relative to `./data`.

### Input Format
Input must contain the description of a nucleic acid sequence of RNA in the following format.

The first and only line must contain a string $s$ $(s_i \in \{A,C,G,U\})$ — the nucleic acid sequence.

### Output Format
The output will contain the description of the maximum matching.

The first line will contain a single integer $m$ — the cardinality of the maximum matching in the sequence. Each of the next $m$ lines will contain two integers — the indices of the base pairs in the matching.

## Visualization
Python script to run the app against some input and plot a graph with [matplotlib](https://matplotlib.org/) using the output.

```sh
cd scripts
./run.py [inputfile]
```

**Note:** `inputfile` may also be relative to `./data`.

## Testing
This project uses [GoogleTest](https://github.com/google/googletest) for its unit tests and [GoogleBenchmark](https://github.com/google/benchmark) for benchmarking.

#### Testing `RNA::Predictor::find_max_matching()` against various cases:
```sh
cd build
ctest -R PredictorTests -j6
```

#### Benchmarking against varying input sizes:
```sh
./bin/bench --benchmark_counters_tabular=true
```
<details>
<summary>Benchmark output</summary>

```python
2022-04-18T03:23:59+05:30
Running ./bin/bench
Run on (12 X 3000 MHz CPU s)
CPU Caches:
  L1 Data 32 KiB (x6)
  L1 Instruction 32 KiB (x6)
  L2 Unified 512 KiB (x6)
  L3 Unified 4096 KiB (x2)
Load Average: 0.91, 0.73, 0.65
***WARNING*** CPU scaling is enabled, the benchmark real time measurements may be noisy and will incur extra overhead.
-------------------------------------------------------------------
Benchmark              Time             CPU   Iterations          n
-------------------------------------------------------------------
BM_Random/4         1369 ns         1323 ns       526810          4
BM_Random/54      150777 ns       150695 ns         4718         54
BM_Random/104    1046964 ns      1046092 ns          672        104
BM_Random/154    3359127 ns      3356692 ns          203        154
BM_Random/204    7698778 ns      7692782 ns           90        204
BM_Random/254   14885152 ns     14874592 ns           47        254
BM_Random/304   25571931 ns     25553385 ns           27        304
BM_Random/354   40506849 ns     40484808 ns           17        354
BM_Random/404   59950695 ns     59919766 ns           11        404
BM_Random/454   85870123 ns     85813275 ns            8        454
--------------------------------------------------------
Benchmark              Time             CPU   Iterations
--------------------------------------------------------
BM_Random_BigO       0.91 N^3        0.91 N^3
BM_Random_RMS          1 %             1 %
```
</details>

</br>

*This page uses math latex formatting. Download the [extension](https://chrome.google.com/webstore/detail/github-math-display/cgolaobglebjonjiblcjagnpmdmlgmda) to render it.*