#include <chrono>
#include <iostream>

#include "MatrixLib.h"

Eigen::MatrixXcd generateTransferMatrix(size_t d, size_t chi) {
    std::vector<Eigen::MatrixXcd> matrixAs = splitMatrix(generatePartialUnitary(d, chi));
    Eigen::MatrixXcd res(chi * chi, chi * chi);
    res.setZero();

    for (auto& i : matrixAs) {
        res += kroneckerProduct(i, i.conjugate());
    }

    return res;
}

int main() {
    /*
    size_t d, chi, ns, it;
    std::cout << "Dimension of each Hilbert space: ";
    std::cin >> d;
    std::cout << "Bond dimension of translation-invariant MPS: ";
    std::cin >> chi;
    std::cout << "Number of sites: ";
    std::cin >> ns;
    std::cout << "Size of ensemble: ";
    std::cin >> it;

    std::complex<double> normSum;

    auto start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < it; i++) {
        Eigen::MatrixXcd transferMat = generateTransferMatrix(d, chi);
        normSum += mat_power(transferMat, ns).trace();
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;

    std::cout << "Average Norm = " << normSum / (double)it << std::endl;
    std::cout << "Calculation time (s) = " << duration.count() << std::endl;*/

    size_t d, chi, it;
    std::cout << "Dimension of each Hilbert space: ";
    std::cin >> d;
    std::cout << "Bond dimension of translation-invariant MPS: ";
    std::cin >> chi;
    std::cout << "Size of ensemble: ";
    std::cin >> it;

    auto start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < it; i++) {
        Eigen::MatrixXcd transferMat = generateTransferMatrix(d, chi);
        Eigen::ComplexEigenSolver<Eigen::MatrixXcd> solver(transferMat);
        Eigen::VectorXcd eigenvalues = solver.eigenvalues();
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;

    std::cout << "Calculation time (s) = " << duration.count() << std::endl;

    return 0;
}