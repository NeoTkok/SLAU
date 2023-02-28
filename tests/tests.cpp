#include <gtest/gtest.h>
#include "Solver_Mat3Dig.hpp"
#include "Class_CSR.hpp"
#include <iostream>


TEST(CLASS_CSR, operator) {
    std::vector<int> col = {0,1,2,0,1,2,0,1,2};
    std::vector<int> row = {0,3,6,9};
    std::vector<double> v = {0,1,2,1,2,3,2,3,4};
    CSR Matrix(col, row, v);
    Matrix *= 10;
    for (int i = 0; i < 3; ++i)
        for(int j = 0; j < 3; ++j)
            ASSERT_DOUBLE_EQ(Matrix(i,j), 10*(i+j));
}

TEST(CLASS_CSR, get){
    std::vector<int> col = {0,1,3,2,1,3};
    std::vector<int> row = {0,3,4,6};
    std::vector<double> v = {1,2,3,4,1,11};
    CSR Matrix(col, row, v);
    for(int z = 0; z < 6; ++z)
        ASSERT_DOUBLE_EQ(Matrix(Matrix.get_i(z),Matrix.get_j(z)), v[z]);
}

TEST(CLASS_CSR, umn_vect){
    std::vector<int> col = {0,1,2,0,1,2,0,1,2};
    std::vector<int> row = {0,3,6,9};
    std::vector<double> v = {0,1,2,1,2,3,2,3,4};
    CSR Matrix(col, row, v);
    std::vector<double> d = {1,2,3};
    std::vector<double> solve = {8, 14, 20};

    for (int i = 0; i < 3; ++i)
        ASSERT_DOUBLE_EQ((Matrix*d)[i], solve[i]);
}



TEST(solver, slae1) {
    int n = 3;
    std::vector<double> a = {1., 1.};
    std::vector<double> b = {4., 4., 4.};
    std::vector<double> c = {1., 1.};
    std::vector<double> d = {28., 7., 56.};
        
    Mat3D A(n,a,b,c);

    auto x = solve(A,d);

    ASSERT_DOUBLE_EQ(x[0], 8.);
    ASSERT_DOUBLE_EQ(x[1], -4.);
    ASSERT_DOUBLE_EQ(x[2], 15.);
}

TEST(solver, slae2) {
    int n = 3;
    std::vector<double> a =  {1., 2.};
    std::vector<double> b =  {11., 12., 13.};
    std::vector<double> c = {1., 2.};
    std::vector<double> d =  {10., 10., 10.};
    
    Mat3D A(n,a,b,c);

    auto x = solve(A,d);
    ASSERT_DOUBLE_EQ(x[0], 470./553.);
    ASSERT_DOUBLE_EQ(x[1], 360./553.);
    ASSERT_DOUBLE_EQ(x[2], 370./553.);
}

TEST(solver, slae3) {
    int n = 10;
    std::vector<double> a = {1., 3., 6., 2., 3., 11., 10., 1., 1.};
    std::vector<double> b = {3.3, 5.2, 6.37, 12.03, 3., 5.6, 95., 102., 2., 1.};
    std::vector<double> c = {2., 1., 2., 5., 1., 1., 50., 30., 0.5};
    std::vector<double> d = {1., 1., 1., 1., 1., 1., 1., 1., 1., 1.};

    Mat3D A(n,a,b,c);
    auto x = solve(A,d);

    double x_ref[10] = {362572824060./1551876085453., 355385766055./3103752170906.,
                                    265300269650./1551876085453., -335582640650./1551876085453.,
                                    1598853453829./3103752170906., -350477627981./3103752170906.,
                                    1349332630563./15518760854530., -186783139020./1551876085453.,
                                    1925442363493./4655628256359., 2730185892866./4655628256359.};
    
    for (std::size_t i = 0; i < 10; ++i) {
        ASSERT_DOUBLE_EQ(x[i], x_ref[i]);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
