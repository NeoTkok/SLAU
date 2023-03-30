#ifndef SLAE_DENSE_HPP
#define SLAE_DENSE_HPP

#include <vector>
#include "alg.hpp"

class Dense{
private:
    int M, N;
    std::vector<double> data;
public:
    Dense(int M, int N, const std::vector<double>& v) : M(M), N(N), data(v) {}
    Dense(int M, int N) : M(M), N(N), data(N*M) {}

    double operator()(int j, int k) const{
        return data[j * N + k];
    }
    
    std::vector<double>  operator*(const std::vector<double>& h) const{
        std::vector<double> newDense(N);
        for(int u = 0; u < M; ++u)
            for(int k = 0; k < N; ++k)
                newDense[u] += data[u * N + k] * h[k];            
        return newDense;
    }

    int get_M() const{
        return M;
    }

    int get_N() const{
        return N;
    }

    const std::vector<double> get_data() const{
        return data;
    }

    void edit(int i, int j, double m){
        data[N*i + j] = m;
    }

    Dense operator-=(const Dense& mat){
        for (int j = 0; j < N*M; ++j)
            data[j] -= mat.data[j];    
        return *this;
    }

    Dense operator+=(const Dense& mat){
        for (int j = 0; j < N*M; ++j)
            data[j] += mat.data[j];        
        return *this;
    }

    Dense operator*=(const double d){
        for (int j = 0; j < N*M; ++j)
            data[j] *= d;        
        return *this;
    }

    Dense operator/=(const double d){
        for (int j = 0; j < N*M; ++j)
            data[j] /= d;        
        return *this;
    }
};

const Dense operator+(const Dense& x1, const Dense& x2){
    return Dense(x1.get_M(),x1.get_N(),x1.get_data()+x2.get_data());
}
const Dense operator-(const Dense& x1, const Dense& x2){
    return Dense(x1.get_M(),x1.get_N(),x1.get_data()-x2.get_data());
}

const Dense operator*(const Dense& x, double l){
    return Dense(x.get_M(),x.get_N(),x.get_data()*l);
}

const Dense operator*(double l, const Dense& x){
    return Dense(x.get_M(),x.get_N(),x.get_data()*l);
}

const Dense operator/(const Dense& x, double l){
    return Dense(x.get_M(),x.get_N(),x.get_data()/l);
}


const Dense operator*(const Dense& x1, const Dense& x2){
    Dense x(x1.get_M(), x2.get_N());
    for(int i = 0; i < x1.get_M(); ++i)
        for(int j = 0; j < x2.get_N(); ++j)
            for (int k = 0; k < x1.get_N(); ++k)
                x.edit(i,j,x1(i,k)*x2(k,j)+x(i,j));
    return x;
}

#endif // SLAE_DENSE_HPP

/*
int main(){
    std::vector<double> x = {1,2,3,4,5,6};
    std::vector<double> y = {1,1,1,1,1,1};
    Dense A(3,2,x);
    Dense B(2,3,y);
    Dense C = A *B;
    for (int i = 0; i < C.get_M(); ++i)
    {
        for(int j = 0; j < C.get_N(); ++j)
            std::cout << C(i,j) << " ";
        std::cout << "\n";
    }
    return 0;

}
*/