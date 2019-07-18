
MESSAGE(STATUS "building mer2seq")

add_executable(mer2seq main.cpp)

target_link_libraries(mer2seq 6mer proc m)
