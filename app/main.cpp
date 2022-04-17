/**
 * @file main.cpp
 * @author the-hyp0cr1t3
 * @brief Main app file
 * @date 2022-04-17
 */
#include <RNA/NucleicAcidSequence.hpp>
#include <RNA/Predictor.hpp>
#include <utils/utils.hpp>

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <chrono>

/**
 * @brief Entry point
 *
 * **Usage** `./app [inputfile]`
 *
 * **Example** `./app sample.txt`
 *
 * @param argc The number of commandline arguments
 * @param argv A list of commandline arguments
 * @return nothing
 */
int main(int argc, char *argv[]) {
    if(argc != 2) {
        std::cerr << "Invalid args; Usage: ./app [inputfile]" << std::endl;
        std::cerr << "Eg: ./app sample.txt" << std::endl;
        exit(EXIT_FAILURE);
    }

    // creating input file stream using argv[1] as fname
    std::ifstream fin;
    std::string fname = argv[1];
    utils::open_file(fin, fname);

    // reading input
    RNA::NASeq seq;
    fin >> seq;

    RNA::Predictor predictor(seq);


    // finding max matchings
    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();

      size_t num_matchings = predictor.find_max_matching();

#ifdef RECOVER_MATCHING
      auto matchings = predictor.recover_matchings();
#endif

    std::chrono::high_resolution_clock::duration total_runtime = std::chrono::high_resolution_clock::now() - start;


    // writing output
    std::cout << std::endl << "  " << num_matchings << std::endl;

#ifdef RECOVER_MATCHING
    for(auto [l, r]: matchings)
        std::cout << "  " << l + 1 << ' ' << r + 1 << std::endl;
#endif

    // time taken
    std::cout << std::endl << "  Took "
              << std::chrono::duration<double, std::milli>(total_runtime).count() << " ms" << std::endl;
}