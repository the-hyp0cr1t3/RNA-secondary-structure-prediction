# Unit testing library
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG release-1.11.0)
FetchContent_MakeAvailable(googletest)
# add GoogleTest

# Benchmarking library
FetchContent_Declare(
  googlebenchmark
  GIT_REPOSITORY https://github.com/google/benchmark.git
  GIT_TAG origin/master)
FetchContent_MakeAvailable(googlebenchmark)
# add GoogleBenchmark

add_subdirectory(benchmark)
add_subdirectory(predictor)

# register a test to check behaviour of app on no commandline args
add_test(NAME NoArgs COMMAND app)
set_tests_properties(NoArgs
  PROPERTIES PASS_REGULAR_EXPRESSION "^Invalid"
)