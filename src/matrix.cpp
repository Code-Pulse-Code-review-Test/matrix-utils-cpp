#include "matrix.h"

#include <cmath>
#include <stdexcept>

Matrix::Matrix(std::size_t rows, std::size_t cols, double fill)
    : rows_(rows), cols_(cols), data_(rows * cols, fill) {}

Matrix Matrix::identity(std::size_t n) {
    Matrix m(n, n);
    for (std::size_t i = 0; i < n; ++i) {
        m.at(i, i) = 1.0;
    }
    return m;
}

double& Matrix::at(std::size_t r, std::size_t c) {
    if (r >= rows_ || c >= cols_) {
        throw std::out_of_range("matrix index out of range");
    }
    return data_[r * cols_ + c];
}

double Matrix::at(std::size_t r, std::size_t c) const {
    if (r >= rows_ || c >= cols_) {
        throw std::out_of_range("matrix index out of range");
    }
    return data_[r * cols_ + c];
}

Matrix Matrix::operator+(const Matrix& other) const {
    if (rows_ != other.rows_ || cols_ != other.cols_) {
        throw std::invalid_argument("matrix sizes do not match");
    }
    Matrix result(rows_, cols_);
    for (std::size_t i = 0; i < data_.size(); ++i) {
        result.data_[i] = data_[i] + other.data_[i];
    }
    return result;
}

Matrix Matrix::operator*(const Matrix& other) const {
    if (cols_ != other.rows_) {
        throw std::invalid_argument("cannot multiply these matrices");
    }
    Matrix result(rows_, other.cols_);
    for (std::size_t i = 0; i < rows_; ++i) {
        for (std::size_t j = 0; j < other.cols_; ++j) {
            double sum = 0.0;
            for (std::size_t k = 0; k < cols_; ++k) {
                sum += at(i, k) * other.at(k, j);
            }
            result.at(i, j) = sum;
        }
    }
    return result;
}

Matrix Matrix::operator*(double scalar) const {
    Matrix result(*this);
    for (double& value : result.data_) {
        value *= scalar;
    }
    return result;
}

bool Matrix::operator==(const Matrix& other) const {
    if (rows_ != other.rows_ || cols_ != other.cols_) {
        return false;
    }
    for (std::size_t i = 0; i < data_.size(); ++i) {
        if (std::fabs(data_[i] - other.data_[i]) > 1e-9) {
            return false;
        }
    }
    return true;
}

Matrix Matrix::transpose() const {
    Matrix result(cols_, rows_);
    for (std::size_t i = 0; i < rows_; ++i) {
        for (std::size_t j = 0; j < cols_; ++j) {
            result.at(j, i) = at(i, j);
        }
    }
    return result;
}

Matrix Matrix::minor(std::size_t skipRow, std::size_t skipCol) const {
    Matrix result(rows_ - 1, cols_ - 1);
    std::size_t r = 0;
    for (std::size_t i = 0; i < rows_; ++i) {
        if (i == skipRow) {
            continue;
        }
        std::size_t c = 0;
        for (std::size_t j = 0; j < cols_; ++j) {
            if (j != skipCol) {
                result.at(r, c++) = at(i, j);
            }
        }
        ++r;
    }
    return result;
}

double Matrix::determinant() const {
    if (rows_ != cols_) {
        throw std::invalid_argument("determinant needs a square matrix");
    }
    if (rows_ == 1) {
        return at(0, 0);
    }
    if (rows_ == 2) {
        return at(0, 0) * at(1, 1) - at(0, 1) * at(1, 0);
    }
    double det = 0.0;
    for (std::size_t j = 0; j < cols_; ++j) {
        const double sign = (j % 2 == 0) ? 1.0 : -1.0;
        det += sign * at(0, j) * minor(0, j).determinant();
    }
    return det;
}

std::ostream& operator<<(std::ostream& out, const Matrix& m) {
    for (std::size_t i = 0; i < m.rows(); ++i) {
        for (std::size_t j = 0; j < m.cols(); ++j) {
            out << m.at(i, j) << (j + 1 < m.cols() ? " " : "");
        }
        out << '\n';
    }
    return out;
}
