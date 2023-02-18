#include <vector>
#include <iostream>

class Mat3D{
private:
    const int n;
    std::vector<double> ThreeVector;    
public:
    Mat3D(int N, const std::vector<double>& a,const std::vector<double>& b,
            const std::vector<double>& c) : n(N), ThreeVector(n*3){
    
    ThreeVector[0] = 0;
    for (int i = 1; i < n; ++i)
        ThreeVector[i] = a[i-1];
    for (int i = 0; i < n; ++i)
        ThreeVector[i + n] = b[i];
    for (int i = 0; i < n-1; ++i)
        ThreeVector[i + 2*n] = c[i];
    ThreeVector[3*n - 1] = 0;
    }
    std::vector<double> get_Mat() const{ // получить вектор состоящий(а,b,c)
        return ThreeVector; 
    }
    int get_n() const{ // получить размер 
        return n;
    }
    //friend std::istream& operator>>(std::istream& is, Mat3D &MatIs);

};

