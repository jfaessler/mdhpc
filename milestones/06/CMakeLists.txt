add_executable(06 06.cpp)
target_link_libraries(06 md_lib)
set(assets scaling.sh)

foreach (asset ${assets})
    add_custom_command(
            TARGET 06 POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy
            ${CMAKE_CURRENT_SOURCE_DIR}/${asset}
            ${CMAKE_CURRENT_BINARY_DIR}/${asset}
            COMMENT "Copying ${asset}"
    )
endforeach ()
