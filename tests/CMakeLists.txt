message(STATUS "Davidson diagonalization tests")

file(GLOB TEST_HEADERS "./*.hpp")

add_executable(davdiag_test main.cpp ${TEST_HEADERS})

target_link_libraries(davdiag_test PUBLIC davdiag)

# Gtest
find_package(GTest REQUIRED)
include(GoogleTest)
target_link_libraries(davdiag_test PUBLIC gtest gtest_main)