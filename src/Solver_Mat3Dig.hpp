#ifndef SLAE_MAT3DIG_HPP
#define SLAE_MAT3DIG_HPP


#include "Class_Mat3Dig.hpp"

std::vector<double> solve(const Mat3D& MatSol, const std::vector<double>& d) {
    std::vector<double> alpha(MatSol.get_n());
    std::vector<double> beta(MatSol.get_n());
        
    alpha[0] = (-MatSol.get_Mat()[2*MatSol.get_n()]/MatSol.get_Mat()[MatSol.get_n()]);
    beta[0] = (d[0]/MatSol.get_Mat()[MatSol.get_n()]);
        
    for (int i = 1; i < MatSol.get_n(); ++i){
        alpha[i] = (-MatSol.get_Mat()[MatSol.get_n()*2+i]/(MatSol.get_Mat()[MatSol.get_n()+i]+MatSol.get_Mat()[i]*alpha[i-1]));
        beta[i] = ((d[i]-MatSol.get_Mat()[i]*beta[i-1])/(MatSol.get_Mat()[MatSol.get_n()+i]+MatSol.get_Mat()[i]*alpha[i-1]));
    } 
    for (int i = MatSol.get_n() - 2; i >= 0; --i)
        beta[i] = beta[i+1] * alpha[i] + beta[i]; 
    
    return beta;

}
/* 
std::istream& operator>>(std::istream& is, Mat3D &MatIs)
{
    //std::cout << "Ввод размера матрицы: ";
    is >> MatIs.n;
    MatIs.ThreeVector.resize(MatIs.n*4);
    MatIs.ThreeVector[0] = 0;
    
    //std:: cout << "Ввод нижней диагонали(" << MatIs.n - 1 << " элементов):";
    for (int i = 1; i < MatIs.n; ++i){
        is >> MatIs.ThreeVector[(i)];
    }
    
    //std:: cout << "Ввод главной диагонали(" << MatIs.n << " элементов):";
    for (int i = 0; i < MatIs.n; ++i){
        is >> MatIs.ThreeVector[(i+MatIs.n)];
    }
    
    //std:: cout << "Ввод верхней диагонали(" << MatIs.n - 1 << " элементов):";
    for (int i = 0; i < MatIs.n-1; ++i){
        is >> MatIs.ThreeVector[(i+2*MatIs.n)];
    }
    MatIs.ThreeVector[(3*MatIs.n-1)] = 0;
    
    //std:: cout << "Свободные члены (" << MatIs.n << " элементов):";
    for (int i = 0; i < MatIs.n; ++i){
        is >> MatIs.ThreeVector[(i+3*MatIs.n)];
    };    
    return is;
}


std::ostream& operator<<(std::ostream& os, const std::vector<double>& ShowV){
    for (auto i: ShowV)
        os << i << std::endl;
    return os;
}

*/

#endif 