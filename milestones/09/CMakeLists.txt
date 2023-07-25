add_executable(09 09.cpp)
target_link_libraries(09 md_lib)

set(assets whisker_small.xyz whisker_large.xyz run.sh)

foreach (asset ${assets})
    add_custom_command(
            TARGET 09 POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy
            ${CMAKE_CURRENT_SOURCE_DIR}/${asset}
            ${CMAKE_CURRENT_BINARY_DIR}/${asset}
            COMMENT "Copying ${asset}"
    )
endforeach ()
