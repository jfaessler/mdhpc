# List of headers
set(MD_HEADERS
        atoms.h
        hello.h
        types.h
        verlet.h
        )

# List of implementation files
set(MD_CPP
        hello.cpp
        verlet.cpp
        atoms.cpp lj_direct_summation.cpp lj_direct_summation.h)

# Create a static library
add_library(md_lib STATIC ${MD_CPP} ${MD_HEADERS})

# Include the current directory as INTERFACE
# This will propagate the include to any target linking against the lib
target_include_directories(md_lib INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

# Link against Eigen (automatically sets up the includes)
target_link_libraries(md_lib PUBLIC Eigen3::Eigen)

# Set up MPI includes and library linking
# This also propagates to further targets
if (MPI_FOUND)
    target_compile_options(md_lib PUBLIC -Wall -Wextra -Wpedantic)
    target_compile_definitions(md_lib INTERFACE USE_MPI)
    target_include_directories(md_lib PUBLIC ${MPI_CXX_INCLUDE_DIRS} SYSTEM)
    target_link_libraries(md_lib PUBLIC MPI::MPI_CXX)
endif ()
