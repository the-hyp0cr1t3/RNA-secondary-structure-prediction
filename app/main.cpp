// #undef RECOVER_MATCHING
#include <RNA/NucleicAcidSequence.hpp>
#include <RNA/Predictor.hpp>
#include <utils/utils.hpp>

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <chrono>

int main(int argc, char *argv[]) {
    if(argc != 2) {
        std::cerr << "Invalid args; Usage ./yee fname" << std::endl;
        std::cerr << "eg: ./yee hello-world.yee" << std::endl;
        exit(EXIT_FAILURE);
    }

    // creating input file stream using argv[1] as fname
    std::ifstream fin;
    std::string fname = argv[1];
    utils::open_file(fin, fname);

    RNA::NASeq seq;
    fin >> seq;

    RNA::Predictor predictor(seq);

    // finding max matchings
    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();

      size_t num_matchings = predictor.find_max_matching();

#ifdef RECOVER_MATCHING
      auto matchings = predictor.recover_matching();
#endif

    std::chrono::high_resolution_clock::duration total_runtime = std::chrono::high_resolution_clock::now() - start;

    std::cout << num_matchings << " matching(s) possible" << std::endl;

#ifdef RECOVER_MATCHING
    for(auto [l, r]: matchings)
        std::cout << l << ' ' << r << std::endl;
#endif

    // time taken
    std::cout << std::endl << "Total runtime: "
              << std::chrono::duration<double, std::milli>(total_runtime).count() << std::endl;
}