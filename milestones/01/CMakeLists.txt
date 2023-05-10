# Defining the executable target
add_executable(01 main.cpp)

# Linking against our MD code (propagates all dependencies)
target_link_libraries(01 PUBLIC md_lib)
