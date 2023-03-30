#ifndef ALG_HPP
#define ALG_HPP

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

std::vector<double> operator/(const std::vector<double>& v, double b){
    std::vector<double> y(v.size());
    for(int i = 0; i < v.size(); ++i){
        y[i] = v[i]/b;
    }
    return y;
}
// норма вектора
double Norma(const std::vector<double>& x) {
    double result = 0.0;
    for (int i = 0; i < x.size();++i) {
        result += x[i]*x[i];
    }
    return sqrt(result);
}
#endif 