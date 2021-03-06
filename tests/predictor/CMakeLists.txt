include(GoogleTest)

# macro to add tests and link with googletest
macro(add_gtest_macro TESTNAME FILES LIBRARIES TEST_WORKING_DIRECTORY)
    add_executable(${TESTNAME} ${FILES})
    target_link_libraries(${TESTNAME} gtest gmock gtest_main ${LIBRARIES})

    #   remove set(CMAKE_CXX_STANDARD xx) if you want target-wise standards
    # target_compile_features(${TESTNAME} PRIVATE cxx_std_xx)

    gtest_discover_tests(${TESTNAME}
      WORKING_DIRECTORY ${TEST_WORKING_DIRECTORY}  # only honoured when run with ctest,
      # see https://developercommunity.visualstudio.com/t/cmake-add-test-working-directory-not-honored-with/427600#T-N1244256
    )
    set_target_properties(${TESTNAME} PROPERTIES FOLDER tests)
endmacro()


# register a test linked with google test
add_gtest_macro(
  predictor_test
  predictor_test.cpp
  predictor
  "${CMAKE_SOURCE_DIR}/data"
)