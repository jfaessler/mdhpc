add_subdirectory(mackay)

add_executable(07 07.cpp)
target_link_libraries(07 md_lib)
add_executable(07_timestep 07_timestep.cpp)
target_link_libraries(07_timestep md_lib)

set(assets cluster_923.xyz cluster_3871.xyz clusters.sh)
set(directories generated_clusters)

foreach (asset ${assets})
    add_custom_command(
            TARGET 07 POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy
            ${CMAKE_CURRENT_SOURCE_DIR}/${asset}
            ${CMAKE_CURRENT_BINARY_DIR}/${asset}
            COMMENT "Copying ${asset}"
    )
endforeach ()
foreach (directory ${directories})
    add_custom_command(
            TARGET 07 POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_directory
            ${CMAKE_CURRENT_SOURCE_DIR}/${directory}
            ${CMAKE_CURRENT_BINARY_DIR}/${directory}
            COMMENT "Copying ${directory}"
    )
endforeach ()
