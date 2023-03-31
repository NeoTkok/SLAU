#include <gtest/gtest.h>
#include "Solver_Mat3Dig.hpp"
#include "Solver_ITR.hpp"
#include "alg.hpp"



TEST(Class_Dense, UmnNaSkal) {
    std::vector<double> a = {1., 3., 6., 0.};
    std::vector<double> b = {3., 9., 18., 0.};

    Dense A(2,2,a);
    A *= 3;

    for (int i = 0; i < 2; ++i)
        for(int j = 0; j < 2; ++j)
            ASSERT_DOUBLE_EQ(A(i,j), b[i*2+j]);
}

TEST(Class_Dense, UmnNaSkalCopy) {
    std::vector<double> a = {1., 3., 6., 0.};
    std::vector<double> b = {3., 9., 18., 0.};
    Dense A(2,2,a);

    A = A * 3;
    
    for (int i = 0; i < 2; ++i)
        for(int j = 0; j < 2; ++j)
            ASSERT_DOUBLE_EQ(A(i,j), b[i*2+j]);
}

TEST(Class_Dense, DELenieNaSkal) {

    std::vector<double> a = {1., 3., 6., 0.};
    std::vector<double> b = {3., 9., 18., 0.};

    Dense A(2,2,b);

    A /= 3;
    
    for (int i = 0; i < 2; ++i)
        for(int j = 0; j < 2; ++j)
            ASSERT_DOUBLE_EQ(A(i,j), a[i*2+j]);
}

TEST(Class_Dense, DELenieNaSkalCopy) {
    std::vector<double> a = {1., 3., 6., 0.};
    std::vector<double> b = {3., 9., 18., 0.};

    Dense A(2,2,b);
    A = A / 3;
    
    for (int i = 0; i < 2; ++i)
        for(int j = 0; j < 2; ++j)
            ASSERT_DOUBLE_EQ(A(i,j), a[i*2+j]);
}

TEST(Class_Dense, sum) {
    std::vector<double> a = {1., 3., 6., 0.};
    std::vector<double> b = {3., 9., 18., 1.};
    std::vector<double> c = {4., 12., 24., 1.};
    Dense A(2,2,a);
    Dense B(2,2,b);
    A += B;
    
    for (int i = 0; i < 2; ++i)
        for(int j = 0; j < 2; ++j)
            ASSERT_DOUBLE_EQ(A(i,j), c[i*2+j]);
}

TEST(Class_Dense, sumCopy) {
    std::vector<double> a = {1., 3., 6., 0.};
    std::vector<double> b = {3., 9., 18., 1.};
    std::vector<double> c = {4., 12., 24., 1.};
    Dense A(2,2,a);
    Dense B(2,2,b);
    Dense C = A + B;
    
    for (int i = 0; i < 2; ++i)
        for(int j = 0; j < 2; ++j)
            ASSERT_DOUBLE_EQ(C(i,j), c[i*2+j]);
}

TEST(Class_Dense, minus) {
    std::vector<double> a = {1., 3., 6., 0.};
    std::vector<double> b = {3., 9., 18., 1.};
    std::vector<double> c = {4., 12., 24., 1.};
    Dense A(2,2,c);
    Dense B(2,2,b);
    A -= B;
    
    for (int i = 0; i < 2; ++i)
        for(int j = 0; j < 2; ++j)
            ASSERT_DOUBLE_EQ(A(i,j), a[i*2+j]);
}

TEST(Class_Dense, minusCopy) {
    std::vector<double> a = {1., 3., 6., 0.};
    std::vector<double> b = {3., 9., 18., 1.};
    std::vector<double> c = {4., 12., 24., 1.};
    Dense A(2,2,c);
    Dense B(2,2,b);
    Dense C = A - B;
    
    for (int i = 0; i < 2; ++i)
        for(int j = 0; j < 2; ++j)
            ASSERT_DOUBLE_EQ(C(i,j), a[i*2+j]);
}

TEST(Class_Dense, umnNaVector) {
    std::vector<double> a = {1., 3., 6., 0.};
    std::vector<double> b = {1., 2.};
    std::vector<double> c = {7., 6.};
    Dense A(2,2,a);
    std::vector<double> C = A*b;
    
    for (int i = 0; i < 2; ++i)
            ASSERT_DOUBLE_EQ(C[i], c[i]);
}

TEST(Class_Dense, umn) {
    std::vector<double> a = {1., 3., 6., 0.};
    std::vector<double> b = {1., 2., 3., 4., 5., 6.};
    std::vector<double> c = {13., 17., 21., 6., 12., 18.};
    Dense A(2,2,a);
    Dense B(2,3,b);
    Dense C = A*B;
    
    for (int i = 0; i < 2; ++i)
        for(int j = 0; j < 2; ++j)
            ASSERT_DOUBLE_EQ(C(i,j), c[i*3+j]);
}



//**************************************************************************************
TEST(CLASS_CSR, UmnNaSkal) {
    std::vector<int> col = {0,1,2,0,1,2,0,1,2};
    std::vector<int> row = {0,3,6,9};
    std::vector<double> v = {0,1,2,1,2,3,2,3,4};
    CSR Matrix(3,3, col, row, v);
    Matrix *= 10;
    for (int i = 0; i < 3; ++i)
        for(int j = 0; j < 3; ++j)
            ASSERT_DOUBLE_EQ(Matrix(i,j), 10*(i+j));
}

TEST(CLASS_CSR, get){
    std::vector<int> col = {0,1,3,2,1,3};
    std::vector<int> row = {0,3,4,6};
    std::vector<double> v = {1,2,3,4,1,11};
    CSR Matrix(3, 4, col, row, v);
    for(int z = 0; z < 6; ++z)
        ASSERT_DOUBLE_EQ(Matrix(Matrix.get_i(z),Matrix.get_j(z)), v[z]);
}

TEST(CLASS_CSR, UmnNaVect){
    std::vector<int> col = {0,1,2,0,1,2,0,1,2};
    std::vector<int> row = {0,3,6,9};
    std::vector<double> v = {0,1,2,1,2,3,2,3,4};
    CSR Matrix(3, 3, col, row, v);
    std::vector<double> d = {1,2,3};
    std::vector<double> solve = {8, 14, 20};

    for (int i = 0; i < 3; ++i)
        ASSERT_DOUBLE_EQ((Matrix*d)[i], solve[i]);
}
//*************************************************************************************
TEST(Mat3Dig, slau1) {
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

TEST(Mat3Dig, slau2) {
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

TEST(Mat3Dig, slau3) {
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
//************************************************************************************
TEST(CSR, MetProstIter) {

    std::vector<int> col = {0,1,2,0,1,2,0,1,2};
    std::vector<int> row = {0,3,6,9};
    std::vector<double> v = {10, -2, 6, 3, 8, -1, 1, 2, 1};
    std::vector<double> b = {1, 2, 3};
    CSR A(3, 3, col, row, v);

    std::vector<double> x0 = {0., 0., 0.};
    std::vector<double> x = {-25/24., 11/12., 53/24.};
    std::vector<double> x_ref = MPI(A, b, x0, 0.01, 1000);
    for (std::size_t i = 0; i < 3; ++i) {
        EXPECT_NEAR(x[i], x_ref[i], 1e-4);
    }
}

TEST(CSR, Yakobi) {
    std::vector<int> col = {0,1,2,0,1,2,0,1,2};
    std::vector<int> row = {0,3,6,9};
    std::vector<double> v = {10, -2, 6, 3, 8, -1, 1, 2, 1};
    std::vector<double> b = {1, 2, 3};
    CSR A(3,3, col, row, v);

    std::vector<double> x0 = {0., 0., 0.};
    std::vector<double> x = {-25/24., 11/12., 53/24.};
    std::vector<double> x_ref = Yakobi(A, b, x0, 1000);
    for (std::size_t i = 0; i < 3; ++i) {
        EXPECT_NEAR(x[i], x_ref[i], 1e-5);
    }
}

TEST(DENSE, MetProstIter) {

    std::vector<double> v = {10, -2, 6, 3, 8, -1, 1, 2, 1};
    std::vector<double> b = {1, 2, 3};
    
    Dense A(3, 3, v);

    std::vector<double> x0 = {0., 0., 0.};
    std::vector<double> x = {-25/24., 11/12., 53/24.};
    std::vector<double> x_ref = MPI(A, b, x0, 0.01, 1000);
    
    for (std::size_t i = 0; i < 3; ++i) {
        EXPECT_NEAR(x[i], x_ref[i], 1e-4);
    }
}

TEST(DENSE, MetGAUSZEY) {

    std::vector<double> v = {10, -2, 6, 3, 8, -1, 1, 2, 1};
    std::vector<double> b = {1, 2, 3};
    
    Dense A(3, 3, v);

    std::vector<double> x0 = {0., 0., 0.};
    std::vector<double> x = {-25/24., 11/12., 53/24.};
    std::vector<double> x_ref = GausZ(A, b, x0, 100);
    
    for (std::size_t i = 0; i < 3; ++i) {
        EXPECT_NEAR(x[i], x_ref[i], 1e-4);
    }
}



int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
