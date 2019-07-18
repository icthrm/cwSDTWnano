
MESSAGE(STATUS "building main")

add_executable(sat-align main.cpp)
target_link_libraries(sat-align proc 6mer m wavelib)
