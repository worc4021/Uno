
add_executable(nlp_test src/hs71.cpp)
target_link_libraries(nlp_test uno)
target_include_directories(nlp_test PRIVATE include)
set_target_properties(nlp_test PROPERTIES 
      POSITION_INDEPENDENT_CODE ON
      MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")