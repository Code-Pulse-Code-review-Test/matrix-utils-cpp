#include <iostream>

#include "matrix.h"

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
    std::cout << "A^-1 =\n" << a.inverse();

    // 4x + 7y = 18, 2x + 6y = 14
    const std::vector<double> x = a.solve({18, 14});
    std::cout << "x = " << x[0] << ", y = " << x[1] << '\n';
    return 0;
}
