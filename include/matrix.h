#pragma once

#include <cstddef>
#include <ostream>
#include <vector>

class Matrix {
public:
    Matrix(std::size_t rows, std::size_t cols, double fill = 0.0);

    static Matrix identity(std::size_t n);

    std::size_t rows() const { return rows_; }
    std::size_t cols() const { return cols_; }

    double& at(std::size_t r, std::size_t c);
    double at(std::size_t r, std::size_t c) const;

    Matrix operator+(const Matrix& other) const;
    Matrix operator*(const Matrix& other) const;
    Matrix operator*(double scalar) const;
    bool operator==(const Matrix& other) const;

    Matrix transpose() const;
    double determinant() const;
    double trace() const;

    // throws std::domain_error when the matrix is singular
    Matrix inverse() const;
    // x such that (*this) * x == b
    std::vector<double> solve(const std::vector<double>& b) const;

private:
    std::size_t rows_;
    std::size_t cols_;
    std::vector<double> data_;

    Matrix minor(std::size_t skipRow, std::size_t skipCol) const;
    void swapRows(std::size_t a, std::size_t b);
};

std::ostream& operator<<(std::ostream& out, const Matrix& m);
