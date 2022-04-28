#include <NucleicAcidSequence.hpp>

bool RNA::matches(RNA::BASE x, RNA::BASE y) {
    switch(x) {
        case RNA::BASE::A :
            return y == RNA::BASE::U;

        case RNA::BASE::C :
            return y == RNA::BASE::G;

        case RNA::BASE::G :
            return y == RNA::BASE::C;

        case RNA::BASE::U :
            return y == RNA::BASE::A;

        default:
            return false;
    }
}


RNA::NASeq::NASeq(const std::string &seq) {
    try {
        RNA::NASeq::validate(seq);
    } catch (const std::runtime_error &e) {
        throw "Invalid nucleic acid base sequence provided";
    }

    build(seq);
}

std::istream &RNA::operator >> (std::istream &is, RNA::NASeq &seq) {
    std::string s;
    is >> s;

    try {
        RNA::NASeq::validate(s);
    } catch (const std::runtime_error &e) {
        throw "Invalid nucleic acid base sequence read";
    }

    seq.build(s);

    return is;
}


void RNA::NASeq::build(const std::string &seq) {
    size = 0;
    for(char c: seq) {
        RNA::BASE b;
        b = c == 'A'? RNA::BASE::A
            : c == 'C'? RNA::BASE::C
            : c == 'G'? RNA::BASE::G
            : RNA::BASE::U;
        sequence[size++] = b;
    }
}

bool RNA::NASeq::validate(const std::string &seq) {

    if(seq.empty())
        throw std::runtime_error("Empty nucleic acid base sequence");

    for(char c: seq)
        if(c != 'A' and c != 'C' and c != 'G' and c != 'U')
            throw std::runtime_error("Invalid nucleic acid base");

    return true;
}

RNA::BASE &RNA::NASeq::operator [] (size_t idx) {
    return sequence[idx];
}

const RNA::BASE &RNA::NASeq::operator [] (size_t idx) const {
    return sequence[idx];
}

size_t RNA::NASeq::length() const {
    return size;
}
