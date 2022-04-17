add_executable(app
  main.cpp
  "${CMAKE_SOURCE_DIR}/src/utils/utils.cpp"

  # headers are optional, but won't show up nicely in IDEs unless listed
  "${CMAKE_SOURCE_DIR}/include"
  "${PROJECT_BINARY_DIR}/include"
)

target_include_directories(app PRIVATE
  "${CMAKE_SOURCE_DIR}/include"
  "${PROJECT_BINARY_DIR}/include"
)

target_link_libraries(app PRIVATE predictor)

set_target_properties(app
  PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin"
)