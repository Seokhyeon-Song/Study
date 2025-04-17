#pragma once

#include <complex>
#include <vector>

using Complex = std::complex<double>;

Complex trace(const std::vector<Complex> &A, const size_t dim);
void matrixPower(size_t n, const std::vector<Complex> &A, std::vector<Complex> &result,
                 const size_t dim);
void complexEigenValues(const std::vector<Complex> &A, std::vector<Complex> &w, const size_t n);

std::vector<Complex> generateTransferMatrix(const size_t d, const size_t chi);
std::vector<Complex> generateTransferMatrix2(const size_t d, const size_t chi);