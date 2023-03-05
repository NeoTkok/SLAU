#include "Class_CSR.hpp"
#include <iostream>
#include <vector>
#include <cmath>

//сложение и вычетание stl-ных векторов
std::vector<double> sum_vect(const std::vector<double>& x1, const std::vector<double>& x2){
    std::vector<double> y(x1.size());
    for(int i = 0; i < x1.size(); ++i){
        y[i] = x1[i] + x2[i];
    }
    return y;
}

//
std::vector<double> minus_vect(const std::vector<double>& x1, const std::vector<double>& x2){
    std::vector<double> y(x1.size());
    for(int i = 0; i < x1.size(); ++i){
        y[i] = x1[i] - x2[i];
    }
    return y;
}

// их же произветедие на скаляр
std::vector<double> compos(const std::vector<double> v, const double a){
    std::vector<double> V = v;
    for (int i = 0; i < V.size(); ++i)
        V[i] *= a;
    return V;
} 

// метод Якоби
std::vector<double> Yakobi(const CSR& MATRIX, const std::vector<double> B, const std::vector<double>& X, int N){
    std::vector<double> interX = X;
    for(int i = 0; i < N; ++i)
        interX = (MATRIX.D()).obrdiag() * minus_vect(B,MATRIX.LplusU()*interX);
    return interX;
}

// Метод простых итараций
std::vector<double> MPI(CSR MATRIX, const std::vector<double> B, const std::vector<double>& X, const double t, int N){
    std::vector<double> interX = X;
    for(int i = 0; i < N; ++i){
        interX = minus_vect(interX, compos(minus_vect(MATRIX*interX,B),t));
    }
    return interX;
}

// норма вектора
double ComputeNorm(const std::vector<double>& x) {
    double result = 0.0;
    for (int i = 0; i < x.size();++i) {
        result += x[i]*x[i];
    }
    return sqrt(result);
}

// функция для 3 задания
int NNN(CSR MATRIX, const std::vector<double> B, const std::vector<double>& X, const double t, const std::vector<double>& RESHENIE){
    std::vector<double> interX = X;
    int i = 0;
    for(;;++i){
        interX = minus_vect(interX, compos(minus_vect(MATRIX*interX,B),t));
     //   std::cout << ComputeNorm(minus_vect(interX, RESHENIE)) << "\n";
        if (ComputeNorm(minus_vect(interX, RESHENIE)) < 1e-12)
            break;
    }
    return i;
}