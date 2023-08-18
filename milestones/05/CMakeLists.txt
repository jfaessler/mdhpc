add_executable(05 05.cpp)
target_link_libraries(05 md_lib)

set(assets scaling.sh)

foreach (asset ${assets})
    add_custom_command(
            TARGET 05 POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy
            ${CMAKE_CURRENT_SOURCE_DIR}/${asset}
            ${CMAKE_CURRENT_BINARY_DIR}/${asset}
            COMMENT "Copying ${asset}"
    )
endforeach ()
