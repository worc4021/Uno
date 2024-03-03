project(nlp_test)

add_executable(nlp_test src/hs71.cpp)
target_link_libraries(nlp_test Uno::uno)
set_property(TARGET nlp_test PROPERTY
            MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
target_include_directories(nlp_test PRIVATE include)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/../uno.options DESTINATION ${CMAKE_CURRENT_BINARY_DIR})