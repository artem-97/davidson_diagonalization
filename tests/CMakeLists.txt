# get_filename_component(LIB_PATH "." ABSOLUTE)
get_filename_component(LIB_INCLUDE_PATH "../src" ABSOLUTE)

# file(GLOB_RECURSE LIB_CXX_SOURCES ${BIN_PATH}/*.cpp)
# file(GLOB_RECURSE LIB_HEADERS ${BIN_PATH}/*.hpp)
add_executable(davdiag_test main.cpp)
target_link_libraries(davdiag_test davdiag)
target_include_directories(davdiag PUBLIC ${LIB_INCLUDE_PATH})
