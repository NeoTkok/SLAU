#ifndef SLAE_GMRES_HPP
#define SLAE_GMRES_HPP


#include <vector>

#include "Class_Dense.hpp"
#include "alg.hpp"



std::pair<Dense,Dense> Arnoldy(const Dense& A, const std::vector<double> B, const std::vector<double>& X, const int n_i){
    std::vector<double> r = A * X - B;
    
    int N = A.get_M();
    Dense H(n_i + 1, n_i);
    std::vector<std::vector<double>> v(N + 1);

    v[0] = r / Norma_2(r);
    for (int j = 0; j < n_i; ++j){
        v[j + 1] = A * v[j];
        for(int k = 0; k <= j; ++k){
            H.edit(k, j, v[j + 1] * v[k]);
            v[j + 1] = v[j + 1] - H(k,j) * v[k];            
        }   
        H.edit(j + 1, j, Norma_2(v[j + 1]));
        v[j + 1] = v[j + 1] / H(j + 1, j); 
    }
    std::vector<double> w(n_i*N);
    for (int i = 0; i < n_i; ++i)
        for (int j = 0; j < N; ++j)
            w[i+(n_i)*j] = (v[i])[j];
    Dense W(N,n_i,w);
    std::pair<Dense,Dense> P(H,W);

/*
    std::cout << "+++++++++++++++++++" << std::endl;
    for(int i = 0; i < v.size(); ++i)
    {
        for (auto j: v[i]){
            std::cout << j << " "; 
        }
        std::cout << std::endl;
    }
    std::cout << "+++++++++++++++++++" << std::endl;*/
    return P;
}

// получение синуса и косинуса при i вращении
std::array<double, 2> cos_sin(const Dense& H, const int i){
    std::array<double,2> T = {0, 0};
    double a = H(i, i);
    double b = H(i + 1, i);
    double r = sqrt(a * a + b * b);
    T[0] = b / r; // sin
    T[1] = -a / r; // cos
    return T;
}

// получение матрицы Гивенса при вращении двух соседних элементов
Dense I_E(const int n, const int j, const double C, const double S){//
    Dense E(n,n);
    for(int i = 0; i < n; ++i)
        E.edit(i,i,1.);
    E.edit(j, j, C);
    E.edit(j, j+1, -S);
    E.edit(j+1, j, S);
    E.edit(j+1, j+1, C);
    
    return E;
}


std::pair<Dense,Dense> Qt_R(const Dense& H){
    Dense H_g = H; 
    int size = H.get_N(); // размер 
    std::vector<double> SIN(size);
    std::vector<double> COS(size);
    for(int k = 0; k < size; ++k){
        for(int j = 0; j < k; ++j){
            double d = H_g(j,k);
            H_g.edit(j,k, COS[j] * d - SIN[j] * H_g(j+1, k));
            H_g.edit(j+1,k, SIN[j] * d + COS[j] * H_g(j+1,k));
        }
        std::array<double, 2> z = cos_sin(H_g, k);
        SIN[k] = z[0]; //
        COS[k] = z[1];
        double d = H_g(k,k);
        H_g.edit(k,k, COS[k] * d - SIN[k] * H_g(k+1, k));
        H_g.edit(k+1,k, SIN[k] * d + COS[k] * H_g(k+1,k));
    }

    Dense E = I_E(size + 1,0,1.,0.);
    for(int i = size - 1; i>=0 ; --i){
        Dense Ei = I_E(size+1, i, COS[i],SIN[i]);
        E = E * Ei;
        
    }
    // таким образом мы получили 2 матрицы Q^T = E  и  R = H_g
    std::pair<Dense,Dense> Z(E,H_g);
    return Z;
}

std::vector<double> gmres(const Dense& A, const std::vector<double>& b,
            const std::vector<double>& x0, const int n_i, const double eps){
    int count = 0;
    std::vector<double> X_i = x0;
    std::vector<double> r = A * x0 - b;

    while (Norma_2(r) > eps){
        std::pair<Dense,Dense> H_W =  Arnoldy(A,b,X_i, n_i);
        Dense H = H_W.first;
        Dense W = H_W.second;
        int N = H.get_N();
        std::pair<Dense,Dense> Z = Qt_R(H);
        //  Z.first = Q^t
        //  Z.second = R

        Dense R_new(n_i, n_i);
        for(int i = 0; i < n_i; ++i)
            for(int j = i; j < n_i; ++j )
                R_new.edit(i,j,Z.second(i,j));
        // получили новую R__i_i 

        Dense Q_new = Z.first;
        // получили новую Q^t__i+1_i+1


        std::vector<double> e1(n_i);
        e1[0] = 1;
        double w = Norma_2(r);
        std::vector<double> z = Z.first * e1 * w;
        std::vector<double> z_new(n_i);
        for(int i = 0; i < n_i; ++i)
            z_new[i] = z[i];
        // получили z__i
    
        // таким образом надо решить обратным методом гаусса Ry=Z
        // т е в моих переменных R_new * y = z_new
        std::vector<double> Yiter(n_i);

        for(int i = n_i - 1; i>=0; --i){
            double sum = 0.;
            for(int j = n_i - 1; j > i; --j)
                sum += Yiter[j] * R_new(i,j);
            Yiter[i] = (z_new[i]-sum) / R_new(i,i);
        }
        X_i = X_i - W * Yiter;
        r = A * X_i - b;
        count++;
    }
    std::cout << "колличество срабатывания проекционного метода = " << count << std::endl;
    return X_i;
}




#endif // SLAE_GMRES_HPP
