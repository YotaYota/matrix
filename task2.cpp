// Author: Johannes Jeppsson

using namespace std;

#include <iostream>
#include "r8lib.h"
#include "r8mat_expm1.h"
#include "matrix.hpp"

Matrix eexp(Matrix A, double tol = 1e-10) {
  int size = A.get_dim();
  double norm;
  norm = 2*tol;
  Matrix term(size);
  Matrix sum(size);
  int n;
  for (n = 1; norm > tol; n++) {
    term *= A/n;
    sum += term;
    norm = term.frobenius_norm();
  }
  cout << " iterations: " << n << endl;
  return sum;
}

int main() {
  Matrix A = Matrix(2);
  cout << " Enter elements (column major): " << endl;
  A.fill_matrix();
  cout << " Matrix A: " << endl;
  A.print_matrix();
  Matrix exp_A = eexp(A, 1e-40);
  cout << " e^A: " << endl;
  exp_A.print_matrix();

  int n = A.get_dim();
  Matrix exp_A_mat = Matrix(n, r8mat_expm1(n, A.get_array()));
  cout << " e^A matlab: " << endl;
  exp_A_mat.print_matrix();

  double norm_A = exp_A.frobenius_norm();
  double norm_A_mat = exp_A_mat.frobenius_norm();
  cout << " norm e^A | norm e^A matlab  <-> ";
  cout << norm_A << " | " << norm_A_mat << endl;
}

