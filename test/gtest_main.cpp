#include <gtest/gtest.h>
#include <iostream>

// the main funciton for gtest
int main(int argc, char* argv[]) {
    int result = 0;

    // set up gtest
    ::testing::InitGoogleTest(&argc, argv);

    // running tests
    result = RUN_ALL_TESTS();

    // return good status no matter what
    return result;
}
