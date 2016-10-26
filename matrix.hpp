// matrix.hpp

#ifndef MATRIX_HPP
#define MATRIX_HPP

// Square matrix
class Matrix {
  private:
    unsigned int n;
    double *A;
  public:
    Matrix(unsigned int);
    Matrix(const Matrix&);
    Matrix(const unsigned int, const double*);

    Matrix& operator=(const Matrix&);
    Matrix& operator+=(const Matrix&);
    Matrix& operator*=(const Matrix&);
    Matrix operator*(const Matrix&);
    Matrix operator/(const unsigned int) const;
    Matrix operator-(const Matrix&) const;

    double frobenius_norm() const;
    void identity();
    void print_matrix() const;
    void fill_matrix();
    double* get_col(const int) const;
    double* get_row(const int) const;
    unsigned int get_dim() const;
    double* get_array() const;
    double dot_product(const int, const double*, const double*);
    void replace_row(const int, const double*);
};

#endif
