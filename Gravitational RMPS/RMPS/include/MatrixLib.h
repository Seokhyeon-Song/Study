#pragma once
#include <complex>
#include <random>
#include <vector>

#define EIGEN_USE_BLAS
#define EIGEN_USE_LAPACKE
#include <Eigen/Dense>

Eigen::MatrixXcd generatePartialUnitary(size_t d, size_t chi) {
    std::mt19937 gen(std::random_device{}());
    std::normal_distribution<> dist(0.0, 1.0);

    // Creating Random Complex Matrix
    Eigen::MatrixXcd Z(d * chi, chi);
    for (size_t i = 0; i < d * chi; ++i)
        for (size_t j = 0; j < chi; ++j)
            Z(i, j) = std::complex<double>(dist(gen), dist(gen));

    // QR Decomposition
    Eigen::HouseholderQR<Eigen::MatrixXcd> qr(Z);
    Eigen::MatrixXcd Q = qr.householderQ() * Eigen::MatrixXcd::Identity(d * chi, chi);

    return Q;
}

Eigen::MatrixXcd kroneckerProduct(const Eigen::MatrixXcd& A, const Eigen::MatrixXcd& B) {
    size_t m = A.rows(), n = A.cols();
    size_t p = B.rows(), q = B.cols();

    Eigen::MatrixXcd result(m * p, n * q);

    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < n; ++j) {
            result.block(i * p, j * q, p, q) = A(i, j) * B;
        }
    }

    return result;
}

std::vector<Eigen::MatrixXcd> splitMatrix(const Eigen::MatrixXcd& A) {
    size_t m = A.cols();
    size_t n = A.rows() / m;
    std::vector<Eigen::MatrixXcd> blocks;
    blocks.reserve(m);

    for (int i = 0; i < m; ++i) {
        blocks.push_back(A.block(i * n, 0, n, n));
    }

    return blocks;
}

Eigen::MatrixXcd mat_power(const Eigen::MatrixXcd& A, size_t n) {
    Eigen::MatrixXcd result = Eigen::MatrixXcd::Identity(A.rows(), A.cols());
    Eigen::MatrixXcd base = A;

    while (n > 0) {
        if (n % 2 == 1)
            result *= base;
        base *= base;
        n /= 2;
    }
    return result;
}