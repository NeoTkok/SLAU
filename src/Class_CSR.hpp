#include <iostream>
#include <vector>


class CSR{
private:
    //const int m;
    //const int n;
    std::vector<int> col = {};
    std::vector<int> row = {} ;
    std::vector<double> val = {};
public:
    CSR(const std::vector<int>& COL, const std::vector<int>& ROW,
        const std::vector<double>& VAL) : val(VAL), col(COL), row(ROW){}
    
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

    std::vector<double> operator*(const std::vector<double>& h){
        std::vector<double> newCSR(row.size());
        for(int u = 0; u < row.size(); ++u)
            for(int k = 0; k < row[u+1]-row[u]; ++k)
                newCSR[u] += val[k+row[u]]*h[col[k+row[u]]];
        return newCSR;
    }
/*
    const CSR& operator*(const CSR& MAT){
        CSR Matrix;
        for ()
            for(int k = 0; k < val.size())
            *this.get_i(k);
            MAT



        return *MAT;
    }

    add_ij(const CSR& MAT,int I; int J, double v){
        for(int k = i;k<;++i)
            if (MAT.get_i(i) == I){
                for(int j = 0; j < J ;++j)


                break;
            }
    }
*/

    const CSR& operator*=(double C){
        for(int i = 0;i < val.size(); ++i)
            val[i] *= C;
        return *this; 
    }

   // CSR operator+(const CSR& Matrix){};

};