#include "RMPSLib.h"

#include <mkl.h>
#include <mkl_lapacke.h>
#include <random>

std::vector<Complex> generateOrthonormalBlock(const size_t d, const size_t chi) {
  const size_t rows = d * chi;
  std::mt19937 gen(std::random_device{}());
  std::normal_distribution<double> dist(0.0, 1.0);
  std::vector<Complex> A(rows * chi);
  for (size_t i = 0; i < rows * chi; ++i) {
    A[i] = Complex(dist(gen), dist(gen));
  }

  std::vector<Complex> tau(chi);
  LAPACKE_zgeqrf(LAPACK_COL_MAJOR, rows, chi, reinterpret_cast<MKL_Complex16 *>(A.data()), rows,
                 reinterpret_cast<MKL_Complex16 *>(tau.data()));
  LAPACKE_zungqr(LAPACK_COL_MAJOR, rows, chi, chi, reinterpret_cast<MKL_Complex16 *>(A.data()),
                 rows, reinterpret_cast<MKL_Complex16 *>(tau.data()));

  std::vector<Complex> blocks(d * chi * chi);
  for (size_t blk = 0; blk < d; ++blk) {
    for (size_t col = 0; col < chi; ++col) {
      for (size_t row = 0; row < chi; ++row) {
        blocks[row + col * chi + blk * chi * chi] = A[(blk * chi + row) + col * rows];
      }
    }
  }

  return blocks;
}

void complexEigenValues(const std::vector<Complex> &A, std::vector<Complex> &w, const size_t n) {
  LAPACKE_zgeev(LAPACK_COL_MAJOR, 'N', 'N', n,
                reinterpret_cast<MKL_Complex16 *>(const_cast<Complex *>(A.data())), n,
                reinterpret_cast<MKL_Complex16 *>(w.data()), nullptr, 1, nullptr, 1);
}

std::vector<Complex> generateTransferMatrix(const size_t d, const size_t chi) {
  const size_t chisq = chi * chi;
  std::vector<Complex> tempMat(chisq * chisq);
  std::vector<Complex> transferMat(chisq * chisq);

  std::vector<Complex> orthoBlock = generateOrthonormalBlock(d, chi);

  cblas_zherk(CblasColMajor, CblasUpper, CblasNoTrans, chisq, d, 1, orthoBlock.data(), chisq, 1,
              tempMat.data(), chisq);
  for (size_t i = 0; i < chisq; ++i) {
    for (size_t j = i + 1; j < chisq; ++j) {
      tempMat[j + i * chisq] = conj(tempMat[i + j * chisq]);
    }
  }

  for (size_t a = 0; a < chi; ++a) {
    for (size_t b = 0; b < chi; ++b) {
      for (size_t c = 0; c < chi; ++c) {
        for (size_t d = 0; d < chi; ++d) {
          // Source indices
          size_t row_S = a * chi + c;
          size_t col_S = b * chi + d;
          size_t idx_S = col_S * chisq + row_S;

          // Target indices
          size_t row_T = a * chi + b;
          size_t col_T = c * chi + d;
          size_t idx_T = col_T * chisq + row_T;

          transferMat[idx_T] = tempMat[idx_S];
        }
      }
    }
  }

  return transferMat;
}