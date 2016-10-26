// matrix.cpp
//
// Class of square matrix with elements ordered according to "column major".
//

#include <cmath>
#include <iostream>
#include "matrix.hpp"

// Constructors

// Constructs the identity matrix of size nxn.
Matrix::Matrix(unsigned int dim) : n(dim) {
  if (n*n > 0) {
    A = new double[n*n];
  }
  identity();
}

Matrix::Matrix(const Matrix& B) : n(B.n) {
  if (n*n > 0) {
    A = new double[n*n];
    for (int i = 0; i < n*n; i++) {
      A[i] = B.A[i];
    }
  }
}

Matrix::Matrix(const unsigned int dim, const double elements[]) : n(dim) {
  A = new double[dim*dim];
  for (int k = 0; k < dim*dim; k++) {
    A[k] = elements[k];
  }
}

// Operators

// Note: If A = B.A, A becomes a pointer to B.A.
Matrix& Matrix::operator=(const Matrix& B) {
  n = B.n;
  for (int k = 0; k < n*n; k++) {
    A[k] = B.A[k];
  }
  return *this;
}

Matrix& Matrix::operator+=(const Matrix& B) {
  if (n == B.n) {
    for (int k = 0; k < n*n; k++) {
      A[k] += B.A[k];
    }
  }
  return *this;
}

// TODO: VERY deep iteration... Use get_* and dot_prod
Matrix Matrix::operator*(const Matrix& B) {
  Matrix C = Matrix(n);
  for (int col = 0; col < n; col++) {
    for (int row = 0; row < n; row++) {
      double dot_prod = 0.0;
      for (int k = 0; k < n; k++) {
        dot_prod += A[row + k*n]*B.A[k + col*n];
      }
      C.A[row + col*n] = dot_prod;
    }
  }
  return C;
}

Matrix& Matrix::operator*=(const Matrix& B) {
  for (int i = 0; i < n; i++) {
    double *row_i_old = get_row(i);
    double *row_i_new = new double[n];
    for (int j = 0; j < n; j++) {
      double *col_b_j = B.get_col(j);
      row_i_new[j] = dot_product(n, row_i_old, col_b_j);
    }
    replace_row(i, row_i_new);
  }
  return *this;
}

Matrix Matrix::operator/(const unsigned int denom) const {
  Matrix B = Matrix(n);
  for (int k = 0; k < n*n; k++) {
    B.A[k] = A[k]/denom;
  }
  return B;
}

Matrix Matrix::operator-(const Matrix& B) const {
  Matrix C = Matrix(n);
  for (int k = 0; k < n*n; k++) {
    C.A[k] = A[k] - B.A[k];
  }
  return C;
}

// Functions

double Matrix::frobenius_norm() const {
  double sum;
  sum = 0.0;
  for (int k = 0; k < n*n; k++) {
    sum += pow(A[k], 2.0);
  }
  sum = pow(sum, 0.5);
  return sum;
}

void Matrix::identity() {
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (i == j) A[i + i*n] = 1;
      else A[i + j*n] = 0;
    }
  }
}

void Matrix::replace_row(const int i, const double new_row[]) {
  for (int k = 0; k < n; k++) {
    A[i + k*n] = new_row[k];
  }
}

double Matrix::dot_product(const int size, const double x[], const double y[]) {
  double sum;
  sum = 0;
  for (int i = 0; i < n; i++) sum += x[i]*y[i];
  return sum;
}

double* Matrix::get_array() const {
  return A;
}

double* Matrix::get_row (const int i) const {
  double *row = new double[n];
  for (int j = 0; j < n; j ++) {
    row[j] = A[i + j*n];
  }
  return row;
}

double* Matrix::get_col (const int i) const {
  double *col = new double[n];
  for (int j = 0; j < n; j ++) {
    col[j] = A[j + i*n];
  }
  return col;
}

unsigned int Matrix::get_dim() const {
  return n;
}

void Matrix::fill_matrix() {
  for (int k = 0; k < n*n; k++) {
    std::cin >> A[k];
  }
}

void Matrix::print_matrix() const {
  for (int i = 0; i < n; i++) {
    std::cout << "[";
    for (int j = 0; j < n; j++) {
      std::cout << " " << A[i + j*n];
    }
    std::cout << "]" << std::endl;
  }
  std::cout << std::endl;
}

