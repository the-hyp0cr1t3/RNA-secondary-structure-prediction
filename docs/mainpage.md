# Documentation for RNA Secondary Structure Prediction {#mainpage}

**Problem:** Given a nucleic acid sequence of RNA, find a [maximum matching](https://en.wikipedia.org/wiki/Matching_(graph_theory)) of \f$ \{A,U\} \f$ or \f$ \{C,G\} \f$ base pairs without knots or sharp turns.

This is a modern C++ implementation that employs (iterative) [dynamic programming](https://en.wikipedia.org/wiki/Dynamic_programming) on intervals to find the cardinality of the maximum matching of base pairs.

The main function is `RNA::Predictor::find_max_matching()`.

**Time Complexity:** \f$ \mathcal{O}(n^3) \f$. <br>
**Space Complexity:** \f$ \mathcal{O}(n^2) \f$.
