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
  std::cout << "Number of sites: ";
  std::cin >> ns;
  std::cout << "Size of ensemble: ";
  std::cin >> it;

  std::ofstream fout("output.txt");

  fout << "Dimension of each Hilbert space: " << d << std::endl;
  fout << "Number of sites: " << ns << std::endl;
  fout << "Size of ensemble: " << it << std::endl;

  auto start = std::chrono::high_resolution_clock::now();

  for (size_t chi = 2; chi <= chiMax; ++chi) {
    size_t chisq = chi * chi;
    std::vector<Complex> transferMat(chisq * chisq);
    std::vector<Complex> powerMat(chisq * chisq);
    double trsum = 0.0, trsqsum = 0.0;
      for (size_t i = 0; i < it; i++) {
        transferMat = generateTransferMatrix2(d, chi);
        matrixPower(ns, transferMat, powerMat, chisq);
        double tr = real(trace(powerMat, chisq));
        trsum += tr;
        trsqsum += tr * tr;
      }
      fout << chi << ", " << (trsqsum/it) - (trsum/it) * (trsum/it) << ", " << (trsum/it) << std::endl;
  }

  auto end = std::chrono::high_resolution_clock::now();
  fout.close();

  std::chrono::duration<double> duration = end - start;
  std::cout << "Calculation time (s) = " << duration.count() << std::endl;

  return 0;
}