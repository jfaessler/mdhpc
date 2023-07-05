add_executable(07 07.cpp)
target_link_libraries(07 md_lib)

set(assets cluster_923.xyz cluster_3871.xyz)

foreach (asset ${assets})
    add_custom_command(
            TARGET 07 POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy
            ${CMAKE_CURRENT_SOURCE_DIR}/${asset}
            ${CMAKE_CURRENT_BINARY_DIR}/${asset}
            COMMENT "Copying ${asset}"
    )
endforeach ()
