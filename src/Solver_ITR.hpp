#ifndef SLAE_ITR_HPP
#define SLAE_ITR_HPP

#include "Class_CSR.hpp"
#include <vector>
#include <alg.hpp>
#include "Class_Dense.hpp"
#include <cmath>

constexpr int N = 64;
// получение корней полинома Чебывышева на [-1;1]
const std::vector<double> chebyshev(){
    std::vector<double> x(N);
    double y = sin(M_PI/N/2);
    x[0] = cos(M_PI/N/2);
    x[N-1] = -x[0];
    double Cos_pi_n = 2 * x[0] * x[0] - 1;
    double Sin_pi_n = 2 * x[0] * y;

    for (int i = 1; i < N/2; ++i){
        x[i] = x[i-1] * Cos_pi_n - y * Sin_pi_n;
        y = y * Cos_pi_n + x[i-1] * Sin_pi_n;
        x[N-i-1] = -x[i];
    }
    return x;
}

// разброс 
constexpr std::array<int, N> razbros(){
    std::array<int,N> A = {0,1};
    for (int i = 1; 2<<i <= N; ++i)
        for (int j = (2<<(i-1))- 1; j >= 0 ; --j){
            A[2*j] = A[j];
            A[2*j+1] = (2<<i)-A[2*j]-1;
        }
    return A;
} 

// ускорение Чебышева для МПИ
std::vector<double> MPI_Cheb(double L_min, double L_max, const CSR& MATRIX,
                const std::vector<double>& B, const std::vector<double>& X){
    const std::vector<double> Chebyshev = chebyshev();
    constexpr std::array<int, N> Razbros = razbros();

    double p = (L_min + L_max)/2; 
    double q = (L_max - L_min)/2;
    
    std::vector<double> interX = X;

    for (int i = 0; i < N; ++i)
        interX = interX - 1/(p+q*Chebyshev[Razbros[i]]) * (MATRIX*interX - B);

    return interX;
}

// метод Якоби
std::vector<double> Yakobi(const CSR& MATRIX, const std::vector<double>& B, const std::vector<double>& X, int n){
    std::vector<double> interX = X;
    CSR LU = MATRIX.LplusU();
    CSR obrD = MATRIX.D().obrdiag();
    for(int i = 0; i < n; ++i)
        interX = obrD * (B - LU*interX);
    return interX;
}

// Методы простых итараций
std::vector<double> MPI(const CSR& MATRIX, const std::vector<double>& B, const std::vector<double>& X, const double t, int n){
    std::vector<double> interX = X;
    for(int i = 0; i < n; ++i){
        interX = interX - (MATRIX*interX - B)*t;
    }
    return interX;
}

std::vector<double> MPI(const Dense& MATRIX, const std::vector<double>& B, const std::vector<double>& X, const double t, int n){
    std::vector<double> interX = X;
    for(int i = 0; i < n; ++i){
        interX = interX - (MATRIX*interX - B)*t;
    }
    return interX;
}

// Гаусс Зейдель*************************************************************************************
std::vector<double> GausZ(const Dense& MATRIX, const std::vector<double>& B, const std::vector<double>& X, int N){
    std::vector<double> interX = X;
    for (int i = 0; i < N; ++i){

        std::vector<double> y = interX;
        for(int k = 0; k < MATRIX.get_M(); ++k){
            if (MATRIX(k,k) != 0)
                {
            std::vector<double> Uk(MATRIX.get_N()); // N = M
            std::vector<double> Lk(MATRIX.get_M());
            for(int j = 0; j < k; ++j)
                Lk[j] = MATRIX(k,j);
            for(int j = k+1; j < Uk.size(); ++j)
                Uk[j] = MATRIX(k,j);
        
            interX[k] = (B[k]-Uk*y - Lk*interX)/MATRIX(k,k);            
                }
        }
    }
    return interX;
}

// градиентный спуск

double f_grad(const Dense& A, const std::vector<double>& b, const std::vector<double>& x0, double eps){
    std::vector<double> x = x0;
    std::vector<double> y = x0;
    std::vector<double> Xpred = x0;
    std::vector<double> r = A * x - b;
    std::vector<double> delta = x0;
    double alpha = (r*r)/(r*(A*r));
    double beta = 0;
    for (int i = 0; Norma(r) > eps; ++i){
        x = y - alpha * (r);
        delta = x - Xpred;
        r = A * x - b;
        beta = (alpha * r * (A*r) - r*r) / (delta*(A*r));
        y = x + beta * delta;
        alpha = (beta*delta*(A*r) + r*r) / (r*(A*r));
        Xpred = x;
    }

    return x*(A*x)/2 - b*x;
}

double f_grad(const CSR& A, const std::vector<double>& b, const std::vector<double>& x0, double eps){
    std::vector<double> x = x0;
    std::vector<double> y = x0;
    std::vector<double> Xpred = x0;
    std::vector<double> r = A * x - b;
    std::vector<double> delta = x0;
    double alpha = (r*r)/(r*(A*r));
    double beta = 0;
    for (int i = 0; Norma(r) > eps; ++i){
        x = y - alpha * (r);
        delta = x - Xpred;
        r = A * x - b;
        beta = (alpha * r * (A*r) - r*r) / (delta*(A*r));
        y = x + beta * delta;
        alpha = (beta*delta*(A*r) + r*r) / (r*(A*r));
        Xpred = x;
    }

    return x*(A*x)/2 - b*x;
}











#endif // SLAE_ITR_HPP