#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>

#include "RMPSLib.h"

void sortByReal(std::vector<std::complex<double>> &w) {
  std::sort(w.begin(), w.end(), [](const std::complex<double> &a, const std::complex<double> &b) {
    return std::real(a) > std::real(b);
  });
}

int main() {
  size_t d, chiMax, it, ns;
  std::cout << "Dimension of each Hilbert space: ";
  std::cin >> d;
  std::cout << "Maximum Bond dimension of translation-invariant MPS: ";
  std::cin >> chiMax;
  std::cout << "Size of ensemble: ";
  std::cin >> it;
  std::cout << "Number of sites: ";
  std::cin >> ns;

  std::ofstream fout("output.txt");

  auto start = std::chrono::high_resolution_clock::now();

  for (size_t chi = 1; chi <= chiMax; ++chi) {
    size_t chisq = chi * chi;
    std::vector<Complex> transferMat(chisq * chisq);
    std::vector<Complex> poweredTransferMat(chisq * chisq);
    double sumTrace = 0, sumTraceSq = 0;

    for (size_t i = 0; i < it; i++) {
      transferMat = generateTransferMatrix(d, chi);
      matrixPower(ns, transferMat, poweredTransferMat, chisq);
      double tr = trace(poweredTransferMat, chisq).real();
      sumTrace += tr;
      sumTraceSq += tr * tr;
    }

    sumTrace /= (double)it;
    sumTraceSq /= (double)it;

    fout << chi << ", " << sumTrace << ", " << sumTraceSq - sumTrace * sumTrace << std::endl;
  }

  auto end = std::chrono::high_resolution_clock::now();
  fout.close();

  std::chrono::duration<double> duration = end - start;
  std::cout << "Calculation time (s) = " << duration.count() << std::endl;

  return 0;
}