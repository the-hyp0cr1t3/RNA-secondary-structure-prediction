add_library(predictor STATIC
  Predictor.cpp
  NucleicAcidSequence.cpp

  # headers are optional, but won't show up nicely in IDEs unless listed
  "${CMAKE_SOURCE_DIR}/include/RNA/NucleicAcidSequence.hpp"
  "${CMAKE_SOURCE_DIR}/include/RNA/Predictor.hpp"
)

target_include_directories(predictor
  PUBLIC "${CMAKE_SOURCE_DIR}/include/RNA"
)