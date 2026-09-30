#include "matrix.h"

#include <cmath>
#include <stdexcept>
#include <utility>

namespace {
// pivots smaller than this count as zero
const double kSingularTolerance = 1e-12;
}

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

double Matrix::trace() const {
    if (rows_ != cols_) {
        throw std::invalid_argument("trace needs a square matrix");
    }
    double sum = 0.0;
    for (std::size_t i = 0; i < rows_; ++i) {
        sum += at(i, i);
    }
    return sum;
}

// Gauss-Jordan elimination on [A | I] with partial pivoting
Matrix Matrix::inverse() const {
    if (rows_ != cols_) {
        throw std::invalid_argument("inverse needs a square matrix");
    }
    const std::size_t n = rows_;
    Matrix a(*this);
    Matrix inv = identity(n);

    for (std::size_t col = 0; col < n; ++col) {
        std::size_t pivot = col;
        for (std::size_t r = col + 1; r < n; ++r) {
            if (std::fabs(a.at(r, col)) > std::fabs(a.at(pivot, col))) {
                pivot = r;
            }
        }
        if (std::fabs(a.at(pivot, col)) < kSingularTolerance) {
            throw std::domain_error("matrix is singular");
        }
        a.swapRows(col, pivot);
        inv.swapRows(col, pivot);

        const double scale = a.at(col, col);
        for (std::size_t c = 0; c < n; ++c) {
            a.at(col, c) /= scale;
            inv.at(col, c) /= scale;
        }

        for (std::size_t r = 0; r < n; ++r) {
            if (r == col) {
                continue;
            }
            const double factor = a.at(r, col);
            for (std::size_t c = 0; c < n; ++c) {
                a.at(r, c) -= factor * a.at(col, c);
                inv.at(r, c) -= factor * inv.at(col, c);
            }
        }
    }
    return inv;
}

// goes through the inverse, fine for the small systems in the assignment
std::vector<double> Matrix::solve(const std::vector<double>& b) const {
    if (b.size() != rows_) {
        throw std::invalid_argument("right-hand side has the wrong size");
    }
    Matrix column(rows_, 1);
    for (std::size_t i = 0; i < rows_; ++i) {
        column.at(i, 0) = b[i];
    }
    const Matrix x = inverse() * column;
    std::vector<double> result(rows_);
    for (std::size_t i = 0; i < rows_; ++i) {
        result[i] = x.at(i, 0);
    }
    return result;
}

void Matrix::swapRows(std::size_t a, std::size_t b) {
    if (a == b) {
        return;
    }
    for (std::size_t c = 0; c < cols_; ++c) {
        std::swap(at(a, c), at(b, c));
    }
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
