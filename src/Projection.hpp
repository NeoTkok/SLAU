#ifndef PROJECTION_HPP
#define PROJECTION_HPP

#include <alg.hpp>
#include "Class_Dense.hpp"
#include "Class_CSR.hpp"

#include <vector>


template<typename Matrix>
std::vector<double> CG(const Matrix& A, const std::vector<double>& B,
                        const std::vector<double>& X, const double eps){
    std::vector<double> x = X;
    std::vector<double> r = A*X-B;
    std::vector<double> Rpred = r;
    std::vector<double> d = r;
    double alpha = (r*r)/(d*(A*d));
    double beta = 0;
    while(Norma(r) > eps)
    for (int i = 0; i < B.size()/* i < A.get_M*/; ++i){
        alpha = (r*r)/(d*(A*d));
        x = x - alpha * d;
        Rpred = r;
        r = A*x - B;
        if (Norma(d) == 0)
            break;
        else{
            beta = (r*r)/(Rpred*Rpred); 
            d = r + beta * d; 
        }
    }
    return x;
}

std::vector<double> BiCG(const Dense& A, const std::vector<double>& B,
                        const std::vector<double>& X, const double eps)
{
    std::vector<double> r = A*X-B;
    std::vector<double> r_ = r;
    std::vector<double> iterX = X;
   
    double ro = r_ * r; 
    std::vector<double> p = r;
    std::vector<double> p_ = r_;

    std::vector<double> z = A * p;
    std::vector<double> z_ = A.T() * p_;
    double q  = ro / (p_ * z);
    iterX = iterX - q * p;
    r  = r - q * z;
    r_ = r_ - q * z_;  
    while (Norma_2(r) > eps){
        double ro_new = r_ * r;
        if (ro_new < eps)
            break;
        double t = ro_new/ro;
        p = r + t * p;
        p_ = r_ + t * p_;
        z = A * p;
        z_ = A.T() * p_;
        q  = ro_new / (p_ * z);
        iterX = iterX - q * p;
        r  = r - q * z;
        r_ = r_ - q * z_;  
        
        ro = ro_new;
    }
    
    return iterX;
}

#endif 