#include "Class_CSR.hpp"
#include <iostream>
#include <vector>
#include <cmath>

//сложение stl-ных векторов
std::vector<double> operator+(const std::vector<double>& a, const std::vector<double> b){
    std::vector<double> y(a.size());
    for(int i = 0; i < a.size(); ++i){
        y[i] = a[i] + b[i];
    }
    return y;
}
//вычетание stl-ных векторов
std::vector<double> operator-(const std::vector<double>& a, const std::vector<double> b){
    std::vector<double> y(a.size());
    for(int i = 0; i < a.size(); ++i){
        y[i] = a[i] - b[i];
    }
    return y;
}

// их же произветедие на скаляр
std::vector<double> operator*(const std::vector<double>& v, double b){
    std::vector<double> y(v.size());
    for(int i = 0; i < v.size(); ++i){
        y[i] = v[i]*b;
    }
    return y;
}
std::vector<double> operator*(double b, const std::vector<double>& v){
    std::vector<double> y(v.size());
    for(int i = 0; i < v.size(); ++i){
        y[i] = v[i]*b;
    }
    return y;
}

// метод Якоби
std::vector<double> Yakobi(const CSR& MATRIX, const std::vector<double> B, const std::vector<double>& X, int N){
    std::vector<double> interX = X;
    for(int i = 0; i < N; ++i)
        interX = (MATRIX.D()).obrdiag() * (B - MATRIX.LplusU()*interX);
    return interX;
}

// Метод простых итараций
std::vector<double> MPI(const CSR& MATRIX, const std::vector<double> B, const std::vector<double>& X, const double t, int N){
    std::vector<double> interX = X;
    for(int i = 0; i < N; ++i){
        interX = interX - (MATRIX*interX - B)*t;
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
int NNN(const CSR& MATRIX, const std::vector<double>& B, const std::vector<double>& X, const double t, const std::vector<double>& RESHENIE){
    std::vector<double> interX = X;
    int i = 0;
    for(;;++i){
        interX = interX - (MATRIX*interX - B)*t;
     //   std::cout << ComputeNorm(interX - RESHENIE) << "\n";
        if (ComputeNorm(interX - RESHENIE) < 1e-12)
            break;
    }
    return i;
}