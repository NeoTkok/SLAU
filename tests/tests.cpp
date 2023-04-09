#ifndef SLAE_TESTs_CPP
#define SLAE_TESTS_CPP

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

    std::vector<double> x0 = {0.12, 310., 430.};
    std::vector<double> x = {-25/24., 11/12., 53/24.};
    std::vector<double> x_ref = MPI(A, b, x0, 0.01, 1e-13);
    for (std::size_t i = 0; i < 3; ++i) {
        EXPECT_NEAR(x[i], x_ref[i], 1e-13);
    }
}

TEST(CSR, MetProstIter_cheb ) {

    std::vector<int> col = {0,1,0,1};
    std::vector<int> row = {0,2,4};
    std::vector<double> v = {9, -3, -3, 9};
    std::vector<double> b = {1, 12};
    CSR A(2, 2, col, row, v);

    std::vector<double> x0 = {1230., 0.2325};
    std::vector<double> x = {5/8., 37/24.};
    std::vector<double> x_ref = MPI_Cheb(6., 12., A, b, x0, 1e-25);
    for (std::size_t i = 0; i < 2; ++i) {
        EXPECT_NEAR(x[i], x_ref[i], 1e-25);
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

TEST(CSR, Gauss_Seidel) {
    std::vector<int> col = {0,1,0,1};
    std::vector<int> row = {0,2,4};
    std::vector<double> v = {9, -3, -3, 9};
    std::vector<double> b = {1, 12};
    CSR A(2, 2, col, row, v);

    std::vector<double> x0 = {1230., 0.2325};
    std::vector<double> x = {5/8., 37/24.};
    std::vector<double> x_ref = Gauss_Seidel(A, b, x0, 1e-13);
    for (std::size_t i = 0; i < 2; ++i) {
        EXPECT_NEAR(x[i], x_ref[i], 1e-13);
    }
}

TEST(CSR, SOR) {
    std::vector<int> col = {0,1,0,1};
    std::vector<int> row = {0,2,4};
    std::vector<double> v = {9, -3, -3, 9};
    std::vector<double> b = {1, 12};
    CSR A(2, 2, col, row, v);

    std::vector<double> x0 = {1230., 0.2325};
    std::vector<double> x = {5/8., 37/24.};
    std::vector<double> x_ref = SOR(A, b, x0, 1e-13, 0.5);
    for (std::size_t i = 0; i < 2; ++i) {
        EXPECT_NEAR(x[i], x_ref[i], 1e-13);
    }
}

TEST(CSR, GausSIM) {
    std::vector<int> col = {0,1,0,1};
    std::vector<int> row = {0,2,4};
    std::vector<double> v = {9, -3, -3, 9};
    std::vector<double> b = {1, 12};
    CSR A(2, 2, col, row, v);

    std::vector<double> x0 = {1230., 0.2325};
    std::vector<double> x = {5/8., 37/24.};
    std::vector<double> x_ref = Gaus_SIM(A, b, x0, 1e-13);
    for (std::size_t i = 0; i < 2; ++i) {
        EXPECT_NEAR(x[i], x_ref[i], 1e-13);
    }
}


TEST(DENSE, MetProstIter) {

    std::vector<double> v = {10, -2, 6, 3, 8, -1, 1, 2, 1};
    std::vector<double> b = {1, 2, 3};
    
    Dense A(3, 3, v);

    std::vector<double> x0 = {0.456, 0.23, 120.};
    std::vector<double> x = {-25/24., 11/12., 53/24.};
    std::vector<double> x_ref = MPI(A, b, x0, 0.01, 1e-13);
    
    for (std::size_t i = 0; i < 3; ++i) {
        EXPECT_NEAR(x[i], x_ref[i], 1e-13);
    }
}

TEST(DENSE, Gauss_Seidel) {

    std::vector<double> v = {10, -2, 6, 3, 8, -1, 1, 2, 1};
    std::vector<double> b = {1, 2, 3};
    
    Dense A(3, 3, v);

    std::vector<double> x0 = {0.2, 230., 0.12};
    std::vector<double> x = {-25/24., 11/12., 53/24.};
    std::vector<double> x_ref = Gauss_Seidel(A, b, x0, 1e-13);
    
    for (std::size_t i = 0; i < 3; ++i) {
        EXPECT_NEAR(x[i], x_ref[i], 1e-13);
    }
}

TEST(DENSE, GausSIM) {

    std::vector<double> v = {10, -2, 6, 3, 8, -1, 1, 2, 1};
    std::vector<double> b = {1, 2, 3};
    
    Dense A(3, 3, v);

    std::vector<double> x0 = {0., 0., 0.};
    std::vector<double> x = {-25/24., 11/12., 53/24.};
    std::vector<double> x_ref = Gaus_SIM(A, b, x0, 1e-15);
    
    for (std::size_t i = 0; i < 3; ++i) {
        EXPECT_NEAR(x[i], x_ref[i], 1e-14);
    }
}

TEST(DENSE, SOR) {

    std::vector<double> v = {10, -2, 6, 3, 8, -1, 1, 2, 1};
    std::vector<double> b = {1, 2, 3};
    
    Dense A(3, 3, v);

    std::vector<double> x0 = {0.2, 320., 0.12};
    std::vector<double> x = {-25/24., 11/12., 53/24.};
    std::vector<double> x_ref = SOR(A, b, x0, 1e-15, 0.5);
    
    for (std::size_t i = 0; i < 3; ++i) {
        EXPECT_NEAR(x[i], x_ref[i], 1e-15);
    }
}

TEST(DENSE, CG) {
    
    std::vector<double> v = {10, 0, 0, 1};
    std::vector<double> b = {1, 2};
    
    Dense A(2, 2, v);

    std::vector<double> x0 = {0.3523, 23534.};
    std::vector<double> x = {0.1, 2};
    std::vector<double> x_ref = CG(A, b, x0,1e-15);
    
    for (std::size_t i = 0; i < 2; ++i) {
        EXPECT_NEAR(x[i], x_ref[i], 1e-15);
    }
}

TEST(CSR, CG) {
    std::vector<int> col = {0,1};
    std::vector<int> row = {0,1,2};
    std::vector<double> v = {10, 1};
    std::vector<double> b = {1, 2};
    CSR A(2,2, col, row, v);

    std::vector<double> x0 = {0.3523, 23534.};
    std::vector<double> x = {0.1, 2};
    std::vector<double> x_ref = CG(A, b, x0,1e-15);
    
    for (std::size_t i = 0; i < 2; ++i) {
        EXPECT_NEAR(x[i], x_ref[i], 1e-15);
    }
}

TEST(DENSE, Grad) {
    std::vector<double> v = {10, 0, 0, 1};
    std::vector<double> b = {1, 2};
    
    Dense A(2, 2, v);

    std::vector<double> x0 = {0.1, 0.2};
    std::vector<double> x = {0.1, 2.};
    std::vector<double> x_ref = f_grad(A, b, x0,1e-24);
    
    for (std::size_t i = 0; i < 2; ++i) {
        EXPECT_NEAR(x[i], x_ref[i], 1e-24);
    }
}

TEST(CSR, Grad) {
    std::vector<int> col = {0,1};
    std::vector<int> row = {0,1,2};
    std::vector<double> v = {10, 1};
    std::vector<double> b = {1, 2};
    CSR A(2,2, col, row, v);

    std::vector<double> x0 = {230.1, 10.};
    std::vector<double> x = {0.1, 2};
    std::vector<double> x_ref = f_grad(A, b, x0,1e-25);

    for (std::size_t i = 0; i < 2; ++i) {
        EXPECT_NEAR(x[i], x_ref[i], 1e-24);
    }
}



// ДЛя тестирования самостоялки....
/*

#include <iostream>
#include <fstream>

TEST(CSR, MetProstIterSAM) {
    int n = 289;
    std::vector<int> col(n + 2*n-2 + 2 * (n-17));
    std::vector<int> row(n+1); //
    std::vector<double> v(n + 2*n-2 + 2 * (n-17)); //
    std::vector<double> B(n);
    double a = 12;
    double b = 28;
    int k = 0;
    for(int i = 0; i < n; ++i){
        B[i] = 1;
        for(int j = 0; j < n; ++j){
            if (i-17 == j){
                v[k] = a; 
                col[k] = j; 
                ++k;
                //std::cout << ".12.";
            }
            if (i-1 == j){
                v[k] = a; 
                col[k] = j; 
                ++k;
                //std::cout << ".12.";
            }
            if (i == j){
                v[k] = 2*b; 
                col[k] = j; 
                ++k;
                //std::cout << ".36.";
            }
            if (i+1 == j){
                v[k] = a;
                col[k] = j; 
                ++k;
                //std::cout << ".12.";
            }
            if (i+17 == j){
                v[k] = a; 
                col[k] = j; 
                ++k;
                //std::cout << ".12.";
            }
            //if (i+17 != j && i+1 != j && i != j && i-1 != j && i-17 != j)
                //std::cout << ".0.";
        }
        row[i+1] = k;
        //std::cout<< std::endl;
    }
    double l_max = (b+2*a*cos(M_PI/(n+1)))+(b+2*a*cos(M_PI/(n+1)));
    double l_min = (b+2*a*cos(M_PI*n/(n+1)))+(b+2*a*cos(M_PI*n/(n+1)));
    std::cout<< l_max << "  " << l_min << std::endl;
    CSR A(n, n, col, row, v);
    
    std::vector<double> x0(n);

    


    std::ofstream out;          // поток для записи
    out.open("hello.txt");      // открываем файл для записи


    std::vector<double> z = Gaus_SIMn(A, B, x0, 1e-15);

    for(int i = 0; i < z.size(); ++i){   
        out  << z[i] << ", ";
    }
    std::cout<<std::endl;

    out.close(); 



}  
 */

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

#endif