#include <cassert>
#include <iostream>
#include <stdexcept>

#include "matrix.h"

static void testIdentity() {
    Matrix m(2, 2, 3.0);
    assert(m * Matrix::identity(2) == m);
}

static void testTranspose() {
    Matrix m(2, 3);
    m.at(0, 2) = 5;
    assert(m.transpose().at(2, 0) == 5);
    assert(m.transpose().rows() == 3);
}

static void testDeterminant() {
    Matrix m(3, 3);
    m.at(0, 0) = 2;
    m.at(1, 1) = 3;
    m.at(2, 2) = 4;
    assert(m.determinant() == 24);
}

static void testBadSizes() {
    bool thrown = false;
    try {
        Matrix(2, 3) + Matrix(3, 2);
    } catch (const std::invalid_argument&) {
        thrown = true;
    }
    assert(thrown);
}

int main() {
    testIdentity();
    testTranspose();
    testDeterminant();
    testBadSizes();
    std::cout << "all tests passed\n";
    return 0;
}
