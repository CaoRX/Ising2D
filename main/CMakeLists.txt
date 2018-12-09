aux_source_directory (. MAIN_SRCS)
add_executable (Ising2Dexec ${MAIN_SRCS})
target_link_libraries(Ising2Dexec Lattice)
target_link_libraries(Ising2Dexec Config)