#include <algorithm>
#include <chrono>
#include <iostream>

#include "RMPSLib.h"

void sortByReal(std::vector<std::complex<double>> &w) {
  std::sort(w.begin(), w.end(), [](const std::complex<double> &a, const std::complex<double> &b) {
    return std::real(a) > std::real(b);
  });
}

int main() {
  size_t d, chi, it;
  std::cout << "Dimension of each Hilbert space: ";
  std::cin >> d;
  std::cout << "Bond dimension of translation-invariant MPS: ";
  std::cin >> chi;
  std::cout << "Size of ensemble: ";
  std::cin >> it;

  std::vector<Complex> eigenValues(chi * chi);

  auto start = std::chrono::high_resolution_clock::now();
  for (size_t i = 0; i < it; i++) {
    complexEigenValues(generateTransferMatrix(d, chi), eigenValues, chi * chi);
  }
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> duration = end - start;

  std::cout << "Calculation time (s) = " << duration.count() << std::endl;
}