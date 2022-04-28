/**
 * @file NucleicAcidSequence.hpp
 * @author the-hyp0cr1t3
 * @brief Describes an enum for RNA bases and a class for Nucleic Acid sequences
 * @date 2022-04-17
 */
#pragma once

#include <iostream>
#include <string>

#define MAX_N 501

namespace RNA {

/// An enum for nucleic acid bases
enum BASE {
    A,  ///< Adenine
    C,  ///< Cytosine
    G,  ///< Guanine
    U   ///< Uracil
};

/**
 * @brief Checks if a pair of nucleic acid bases match
 *
 * @param x One nucleic acid base
 * @param y Another nucleic acid base
 * @return `true` if the nucleic acid bases match
 * @return `false` otherwise
 */
bool matches(RNA::BASE x, RNA::BASE y);

/**
 * @brief A class to manage nucleic acid sequences
 */
class NASeq {
    size_t size {0};              ///< The number of nucleic acid bases
    RNA::BASE sequence[MAX_N];    ///< An array storing the sequence of nucleic acid bases

public:

    /**
     * @brief Default constructor
     */
    NASeq() = default;

    /**
     * @brief String constructor
     *
     * Calls `validate()` and then `build()` on \p `seq`
     *
     * @param seq The sequence of nucleic acid bases
     */
    NASeq(const std::string &seq);

    /**
     * @brief Validates if \p `seq` is a proper nucleic acid base sequence
     *
     * @param seq The nucleic acid base sequence
     * @return `true` if \p `seq` is a proper nucleic acid base sequence
     * @return `false` otherwise
     *
     * @throws runtime_error if the sequence is empty
     * @throws runtime_error if the sequence contains any extra symbols (i.e. symbols other than A, C, G, U).
     */
    static bool validate(const std::string &seq);

    /**
     * @brief Builds `sequence` from the nucleic acid base sequence string
     * @pre `validate()` must be called first
     *
     * @param seq The nucleic acid base sequence
     */
    void build(const std::string &seq);

    /**
     * @brief Overloading the [] operator for NASeq
     *
     * @param idx The index of the nucleic acid base
     * @return RNA::BASE& A nucleic acid base
     */
    RNA::BASE &operator [] (size_t idx);

    /**
     * @brief Overloading the [] operator for NASeq
     *
     * @param idx The index of the nucleic acid base
     * @return const RNA::BASE& A nucleic acid base
     */
    const RNA::BASE &operator [] (size_t idx) const;

    /**
     * @brief Returns the length of the sequence
     *
     * @return size_t The length of the sequence
     */
    size_t length() const;

};

/**
 * @brief Overloading the istream >> operator for NASeq
 *
 * @param is A reference to the istream
 * @param s The nucleic acid sequence class
 * @return std::istream& A reference to the istream after the process
 */
std::istream &operator >> (std::istream &is, NASeq &s);

} // RNA