#include <iostream>
#include <vector>


class CSR{
private:
    const int m;
    const int n;
    std::vector<int> col = {};
    std::vector<int> row = {} ;
    std::vector<double> val = {};
public:
    CSR(int M, int N, const std::vector<int>& COL, const std::vector<int>& ROW,
        const std::vector<double>& VAL) : m(M), n(N), val(VAL), col(COL), row(ROW){}

    double operator()(int i, int j) const{
        for (int k = row[i]; k < row[i+1]; ++k)
            if (col[k] == j)
                return val[k]; 
        return 0;
    }

    int get_j (int k) const{
        return col[k];
    }

    int get_i (int k) const{
        int i = 0;
        for (;; ++i)
            if(k < row[i])
                break;
        return i-1;
    }

    const std::vector<double> operator*(const std::vector<double>& h) const{
        std::vector<double> newCSR(m);
        for(int u = 0; u < m; ++u)
            for(int k = 0; k < row[u+1]-row[u]; ++k)
                newCSR[u] += val[k+row[u]]*h[col[k+row[u]]];
        return newCSR;
    }

    CSR L() const{
        CSR l = *this;
        for(int k = 0; k < val.size(); ++k)
            if (get_i(k) <= get_j(k)) 
                l.val[k] = 0;
        return l;
    }
    CSR D() const{
        CSR l = *this;
        for(int k = 0; k < val.size(); ++k)
            if (get_i(k) != get_j(k)) 
                l.val[k] = 0;
        return l;
    }
    CSR U() const{
        CSR l = *this;
        for(int k = 0; k < val.size(); ++k)
            if (get_i(k) >= get_j(k)) 
                l.val[k] = 0;
        return l;
    };

    CSR LplusU() const {
    CSR l = *this;
        for(int k = 0; k < val.size(); ++k)
            if (get_i(k) == get_j(k)) 
                l.val[k] = 0;
        return l;
    }

    CSR obrdiag() const{
    CSR l = *this;
        for(int k = 0; k < val.size(); ++k)
            if(l.val[k] != 0)
                l.val[k] = 1/(l.val[k]);
        return l;
    }
    

    const CSR& operator*=(double C){
        for(int i = 0;i < val.size(); ++i)
            val[i] *= C;
        return *this; 
    }
};