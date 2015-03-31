#include <gtest/gtest.h>
#include <matrix.hpp>

TEST(MatrixTest, MatrixInt) {
    // create small integer matrix
    matrix<int> mat(5,5);

    // check dimensions
    ASSERT_EQ(5, mat.ncols());
    ASSERT_EQ(5, mat.nrows());

    // set some matrix elements:
    mat(0,0) = 1;
    mat(1,1) = 2;
    mat(2,2) = 3;
    mat(3,3) = 4;
    mat(4,4) = 5;

    // check if elements return the same value
    for (int i = 0; i < 5; ++i)
    {
        for (int j = 0; j < 5; ++j)
        {
            if (i == j)
                EXPECT_EQ(i+1, mat(i,j));
            else
                EXPECT_EQ(0, mat(i,j));
        }
    }
}

// TODO: add more tests (different constructors etc)


TEST(MatrixTest, BorderMatrixInt) {
    // create small integer matrix
    int b = 3;
    border_matrix<int> mat(5,5,b);

    // check dimensions
    ASSERT_EQ(5, mat.ncols());
    ASSERT_EQ(5, mat.nrows());
    ASSERT_EQ(b, mat.bordersize());

    // check dims of underlying matrix
    ASSERT_EQ(5+2*b, mat.mat().ncols());
    ASSERT_EQ(5+2*b, mat.mat().nrows());

    // set some matrix elements:
    mat(0,0) = 1;
    mat(1,1) = 2;
    mat(2,2) = 3;
    mat(3,3) = 4;
    mat(4,4) = 5;

    // check if elements return the same value
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < 5; ++j) {
            // check with regular access
            if (i == j)
                EXPECT_EQ(i+1, mat(i,j));
            else
                EXPECT_EQ(0, mat(i,j));

            // check with .at()
            if (i == j)
                EXPECT_EQ(i+1, mat.at(i,j));
            else
                EXPECT_EQ(0, mat.at(i,j));
            // check with raw matrix (offseting by border size
            if (i == j)
                EXPECT_EQ(i+1, mat.mat().at(i+b,j+b));
            else
                EXPECT_EQ(0, mat.mat().at(i+b,j+b));
        }
    }
}
