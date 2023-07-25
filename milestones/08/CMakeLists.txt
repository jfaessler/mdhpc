add_executable(08 08.cpp)
target_link_libraries(08 md_lib)

set(assets cluster_923.xyz cluster_3871.xyz)

foreach (asset ${assets})
    add_custom_command(
            TARGET 08 POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy
            ${CMAKE_CURRENT_SOURCE_DIR}/${asset}
            ${CMAKE_CURRENT_BINARY_DIR}/${asset}
            COMMENT "Copying ${asset}"
    )
endforeach ()
