cmake_minimum_required(VERSION 2.6)

# project settings
project(battle-sim-tests)

FILE(GLOB TEST_FILES *.cpp)
add_executable(test-all ${TEST_FILES} gtest_main.cpp)
target_link_libraries(test-all gtest)
