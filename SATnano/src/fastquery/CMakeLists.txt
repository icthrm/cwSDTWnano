
MESSAGE(STATUS "building main")

add_executable(sat-query main.cpp)
target_link_libraries(sat-query proc 6mer m wavelib)
