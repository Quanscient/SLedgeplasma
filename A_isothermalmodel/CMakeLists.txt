custom_add_executable_from_dir(A_isothermalmodel1D ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(A_isothermalmodel1D sparselizard)
custom_symlink_file(A_isothermalmodel1D ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} "*.msh;*.geo;*.txt;*.csv;*.nas;*.slz")

