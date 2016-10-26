#include <cmath>
#include <iostream>

using namespace std;

// Calculates Euluer's number to the power of x with accuracy at least tol.
// The error is approximated by the N:th term.
double eexp(double x, double tol = 1e-10) {
  double term;
  double sum;
  term = 1.0;
  sum = 1.0;
  for (int n = 1; fabs(term) > tol; n++) {
    term *= x/n;
    sum += term;
  }
  return sum;
}


int main() {
  double x;
  cout << " Enter x: " << flush;
  cin >> x;
  double est = eexp(x);
  double e_pow = exp(x);
  cout << " Approxiamted value: " << est << endl;
  cout << " Built-in value: " << e_pow << endl;
  cout << " Difference: " << fabs(est - e_pow) << endl;
  return 0;
}

