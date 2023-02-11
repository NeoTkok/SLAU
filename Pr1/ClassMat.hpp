#include <vector>
#include <array>
#include <cmath>
#include <stdexcept>
#include <iostream>

class Mat3D{
protected:
    std::vector<float> A;
    int n;
public:
    Mat3D(int N){
        this->n = N;
        A.push_back(0);
        std:: cout << "Ввод нижней диагонали(" << N - 1 << " элементов):";
        for (int i = 1; i < N; ++i){
            float x; std::cin >> x; A.push_back(x);
        }
        std:: cout << "Ввод главной диагонали(" << N << " элементов):";
        for (int i = 0; i < N; ++i){
            float x; std::cin >> x; A.push_back(x);
        }
        std:: cout << "Ввод верхней диагонали(" << N - 1 << " элементов):";
        for (int i = 0; i < N-1; ++i){
            float x; std::cin >> x; A.push_back(x);
        }
        A.push_back(0);
        std:: cout << "Свободные члены (" << N << " элементов):";
        for (int i = 0; i < N; ++i){
            float x; std::cin >> x; A.push_back(x);
        }
    }

    void Reshat(){
        std::vector<double> alpha;
        std::vector<double> beta;
        alpha.push_back(-A[2*n]/A[n]);
        beta.push_back(A[3*n]/A[n]);
        for (int i = 1; i < this->n; ++i){
            alpha.push_back(-A[n*2+i]/(A[n+i]+A[i]*alpha[i-1]));
            beta.push_back((A[3*n+i]-A[i]*beta[i-1])/(A[n+i]+A[i]*alpha[i-1]));
        } 
        for (int i = n - 2; i >= 0; --i)
            beta[i] = beta[i+1] * alpha[i] + beta[i]; 
        for (int i = 0; i < n; ++i)
            std::cout << "x" << i << " = " << beta[i] << std::endl;
    }
};