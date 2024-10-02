custom_add_executable_from_dir(B_selfconsistenttemperature ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(B_selfconsistenttemperature sparselizard)
custom_symlink_file(B_selfconsistenttemperature ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} "*.msh;*.geo;*.txt;*.csv;*.nas;*.slz")

