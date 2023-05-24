include(FetchContent)
include(GoogleTest)

# Fetching Googletest library
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/release-1.12.1.tar.gz
  URL_HASH SHA256=81964fe578e9bd7c94dfdb09c8e4d6e6759e19967e397dbea48d1c10e45d0df2
)

FetchContent_MakeAvailable(googletest)

# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Define headers and test implementation files
set(TESTS_HEADERS
  # empty so far
)

set(TESTS_CPP
  test_hello_world.cpp
  test_verlet.cpp
        test_lj_direct_summation.cpp)

# Creating the test executable
add_executable(tests ${TESTS_CPP} ${TESTS_HEADERS})

# Linking against GTest and our MD code
target_link_libraries(tests PUBLIC md_lib gtest gtest_main)

# Discover Google tests
gtest_discover_tests(tests)

# For tests that do not use GTest
# add_test(NAME <test name> COMMAND <test executable>)
