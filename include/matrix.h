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

private:
    std::size_t rows_;
    std::size_t cols_;
    std::vector<double> data_;

    Matrix minor(std::size_t skipRow, std::size_t skipCol) const;
};

std::ostream& operator<<(std::ostream& out, const Matrix& m);
