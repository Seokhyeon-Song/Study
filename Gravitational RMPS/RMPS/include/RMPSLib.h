#pragma once

#include <complex>
#include <vector>

using Complex = std::complex<double>;

Complex trace(const std::vector<Complex> &A, const size_t dim);
void matrixPower(size_t n, const std::vector<Complex> &A, std::vector<Complex> &result,
                 const size_t dim);
void complexEigenValues(const std::vector<Complex> &A, std::vector<Complex> &w, const size_t n);
void matrixExp(const std::vector<Complex> &A, std::vector<Complex> &expA, size_t chi, const Complex mult);

std::vector<Complex> generateTransferMatrix(const size_t d, const size_t chi);
std::vector<Complex> generateTransferMatrix2(const size_t d, const size_t chi);

void constructReducedDensityMat(const std::vector<Complex> &transferMat,
                                std::vector<Complex> &reducedDensityMat, const size_t np,
                                const size_t nq, const size_t dimsqrt);
double vonNeumannEntropy(const std::vector<Complex> &densityMat, const size_t dim);