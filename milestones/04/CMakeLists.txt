add_executable(04 04.cpp)
target_link_libraries(04 md_lib)
set(assets lj54.xyz)

foreach (asset ${assets})
    add_custom_command(
            TARGET 04 POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy
            ${CMAKE_CURRENT_SOURCE_DIR}/${asset}
            ${CMAKE_CURRENT_BINARY_DIR}/${asset}
            COMMENT "Copying ${asset}"
    )
endforeach ()