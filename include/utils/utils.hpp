#pragma once

#include <fstream>
#include <iostream>

#define EXAMPLES_PATH "data/"

namespace utils {

void open_file(std::ifstream &inFile, const std::string &fname);

} // utils