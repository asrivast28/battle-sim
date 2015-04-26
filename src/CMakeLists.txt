cmake_minimum_required(VERSION 2.6)

# project settings
project(battle-sim)

# TODO: add executables
add_executable(battle main.cpp)


add_executable(test-dmat test_distributed.cpp)
target_link_libraries(test-dmat ${EXTRA_LIBS})
