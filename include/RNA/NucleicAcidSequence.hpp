#pragma once

#include <iostream>
#include <string>
#include <vector>

namespace RNA {

enum BASE { A, C, G, U };

bool matches(RNA::BASE x, RNA::BASE y);

class NASeq {
    std::vector<RNA::BASE> sequence;

public:

    NASeq() = default;

    NASeq(const std::string &seq);

    static bool validate(const std::string &seq);

    void build(const std::string &seq);

    RNA::BASE &operator [] (size_t idx);

    const RNA::BASE &operator [] (size_t idx) const;

    size_t length() const;

};

std::istream &operator >> (std::istream &is, NASeq &s);

} // RNA