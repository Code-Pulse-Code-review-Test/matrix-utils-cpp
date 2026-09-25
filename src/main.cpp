#include <iostream>

#include "matrix.h"

//main
int main() {
    Matrix a(2, 2);
    a.at(0, 0) = 4;
    a.at(0, 1) = 7;
    a.at(1, 0) = 2;
    a.at(1, 1) = 6;

    std::cout << "A =\n" << a;
    std::cout << "A^T =\n" << a.transpose();
    std::cout << "A * I =\n" << a * Matrix::identity(2);
    std::cout << "det(A) = " << a.determinant() << '\n';
    return 0;
}
