#include <vector>
#include <iostream>

class Mat3D{
private:
    int n;
    std::vector<double> ThreeVector;    
public:
    Mat3D(int N, std::vector<double> aa, std::vector<double> bb,
            std::vector<double> cc, std::vector<double> dd) : n(N), ThreeVector(n*4){
    ThreeVector[0] = 0;
    for (int i = 1; i < n; ++i)
        ThreeVector[i]=aa[i-1];
    for (int i = 0; i < n; ++i)
        ThreeVector[i+n]=bb[i];
    for (int i = 0; i < n-1; ++i)
        ThreeVector[i+2*n]=cc[i];
    ThreeVector[3*n-1] = 0;
    for (int i = 0; i < n; ++i)
        ThreeVector[i+3*n]=dd[i];    
    };
    std::vector<double> get_Mat() const{ // получить вектор состоящий(а,b,c,d)
        return ThreeVector; 
    }
    int get_n() const{ // получить размер 
        return n;
    }
    //friend std::istream& operator>>(std::istream& is, Mat3D &MatIs);

};

std::vector<double> solve(const Mat3D& MatSol){
        std::vector<double> alpha(MatSol.get_n());
        std::vector<double> beta(MatSol.get_n());
        
        alpha[0] = (-MatSol.get_Mat()[2*MatSol.get_n()]/MatSol.get_Mat()[MatSol.get_n()]);
        beta[0] = (MatSol.get_Mat()[3*MatSol.get_n()]/MatSol.get_Mat()[MatSol.get_n()]);
        
        for (int i = 1; i < MatSol.get_n(); ++i){
            alpha[i] = (-MatSol.get_Mat()[MatSol.get_n()*2+i]/(MatSol.get_Mat()[MatSol.get_n()+i]+MatSol.get_Mat()[i]*alpha[i-1]));
            beta[i] = ((MatSol.get_Mat()[3*MatSol.get_n()+i]-MatSol.get_Mat()[i]*beta[i-1])/(MatSol.get_Mat()[MatSol.get_n()+i]+MatSol.get_Mat()[i]*alpha[i-1]));
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