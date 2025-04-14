#pragma once

#include <complex>
#include <vector>

using Complex = std::complex<double>;

std::vector<Complex> generateTransferMatrix(const size_t d, const size_t chi);

void complexEigenValues(const std::vector<Complex> &A, std::vector<Complex> &w, const size_t n);