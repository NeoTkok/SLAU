#include "ClassMat.hpp"

int main() {
    int N;
    std::cout << "Ввод размера матрицы: ";
    std::cin >> N;
/*--------------------------------------------------------------------------*/
    Mat3D A(N);
    A.Reshat();
    return 0;
}