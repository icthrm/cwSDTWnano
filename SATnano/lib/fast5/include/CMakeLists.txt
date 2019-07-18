

aux_source_directory(. FAST5_SRCS)

MESSAGE(STATUS "building fast5")

add_library (fast5 STATIC fast5.hpp ${FAST5_SRCS})
target_link_libraries(fast5 hdf5)

set_target_properties(fast5 PROPERTIES LINKER_LANGUAGE Cxx)
set_target_properties(fast5 PROPERTIES OUTPUT_NAME "fast5")
set_target_properties(fast5 PROPERTIES CLEAN_DIRECT_OUTPUT 1)
