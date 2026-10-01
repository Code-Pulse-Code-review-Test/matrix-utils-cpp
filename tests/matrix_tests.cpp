#include <cassert>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

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

static bool approxEqual(double a, double b) {
    return std::fabs(a - b) < 1e-9;
}

static void testTrace() {
    assert(Matrix::identity(3).trace() == 3);
}

static void testInverse() {
    Matrix m(2, 2);
    m.at(0, 0) = 4;
    m.at(0, 1) = 7;
    m.at(1, 0) = 2;
    m.at(1, 1) = 6;
    assert(m * m.inverse() == Matrix::identity(2));
    assert(approxEqual(m.inverse().at(0, 0), 0.6));
}

static void testInverseNeedsRowSwap() {
    // zero in the top left corner
    Matrix m(3, 3);
    m.at(0, 1) = 1;
    m.at(1, 0) = 1;
    m.at(2, 2) = 2;
    assert(m * m.inverse() == Matrix::identity(3));
}

static void testSingular() {
    Matrix m(2, 2, 1.0);
    bool thrown = false;
    try {
        (void)m.inverse();
    } catch (const std::domain_error&) {
        thrown = true;
    }
    assert(thrown);
}

static void testSolve() {
    // 2x + y = 5, x + 3y = 10
    Matrix m(2, 2);
    m.at(0, 0) = 2;
    m.at(0, 1) = 1;
    m.at(1, 0) = 1;
    m.at(1, 1) = 3;
    const std::vector<double> x = m.solve({5, 10});
    assert(approxEqual(x[0], 1));
    assert(approxEqual(x[1], 3));
}

static void testSubtract() {
    Matrix a(2, 2, 5.0);
    Matrix b(2, 2, 2.0);
    assert((a - b) == Matrix(2, 2, 3.0));
    assert((a - a) == Matrix(2, 2));
}

static void testPower() {
    // fibonacci matrix, [[1,1],[1,0]]^10 has F(11) in the corner
    Matrix f(2, 2, 1.0);
    f.at(1, 1) = 0;
    assert(f.power(10).at(0, 0) == 89);
    assert(f.power(0) == Matrix::identity(2));
    assert(f.power(3) == f * f * f);
}

int main() {
    testIdentity();
    testTranspose();
    testDeterminant();
    testBadSizes();
    testTrace();
    testInverse();
    testInverseNeedsRowSwap();
    testSingular();
    testSolve();
    testSubtract();
    testPower();
    std::cout << "all tests passed\n";
    return 0;
}
