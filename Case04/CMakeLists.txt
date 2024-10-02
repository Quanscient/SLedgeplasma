custom_add_executable_from_dir(Case04 ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(Case04 sparselizard)
custom_symlink_file(Case04 ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} "*.msh;*.geo;*.txt;*.csv;*.nas;*.slz")

