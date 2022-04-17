# Documentation for RNA Secondary Structure Prediction {#mainpage}

**Problem:** Given a nucleic acid base RNA sequence, find a maximal matching of \f$ \{A,U\} \f$ or \f$ \{C,G\} \f$ base pairs without knots or sharp turns.

This is a modern C++ implementation that employs (iterative) dynamic programming on intervals <br>
to find the cardinality of the maximum matching of base pairs in \f$ \mathcal{O}(n^3) \f$.

The main function is `RNA::Predictor::find_max_matching()`.