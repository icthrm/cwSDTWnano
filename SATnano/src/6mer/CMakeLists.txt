
aux_source_directory(. MER_SRCS)

MESSAGE(STATUS "building 6mer")

add_library (6mer STATIC ${MER_SRCS})

target_link_libraries(6mer)

set_target_properties(6mer PROPERTIES OUTPUT_NAME "6mer")
set_target_properties(6mer PROPERTIES CLEAN_DIRECT_OUTPUT 1)
