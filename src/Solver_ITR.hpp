#ifndef SLAE_ITR_HPP
#define SLAE_ITR_HPP


#include "Class_CSR.hpp"
#include <vector>
#include <alg.hpp>
#include "Class_Dense.hpp"

// метод Якоби
std::vector<double> Yakobi(const CSR& MATRIX, const std::vector<double>& B, const std::vector<double>& X, int N){
    std::vector<double> interX = X;
    for(int i = 0; i < N; ++i)
        interX = (MATRIX.D()).obrdiag() * (B - MATRIX.LplusU()*interX);
    return interX;
}



// Метод простых итараций
std::vector<double> MPI(const CSR& MATRIX, const std::vector<double>& B, const std::vector<double>& X, const double t, int N){
    std::vector<double> interX = X;
    for(int i = 0; i < N; ++i){
        interX = interX - (MATRIX*interX - B)*t;
    }
    return interX;
}

std::vector<double> MPI(const Dense& MATRIX, const std::vector<double>& B, const std::vector<double>& X, const double t, int N){
    std::vector<double> interX = X;
    for(int i = 0; i < N; ++i){
        interX = interX - (MATRIX*interX - B)*t;
    }
    return interX;
}
//*************************************************************************************
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





// функция для 3 задания
/*
int NNN(const CSR& MATRIX, const std::vector<double>& B, const std::vector<double>& X, const double t, const std::vector<double>& RESHENIE){
    std::vector<double> interX = X;
    int i = 0;
    for(;;++i){
        interX = interX - (MATRIX*interX - B)*t;
     //   std::cout << Norma(interX - RESHENIE) << "\n";
        if (Norma(interX - RESHENIE) < 1e-12)
            break;
    }
    return i;
}
*/


#endif // SLAE_ITR_HPP