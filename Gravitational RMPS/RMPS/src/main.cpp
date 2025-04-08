#include <chrono>
#include <random>

#include "HermitianMatrix.h"

int main() {
    size_t n;
    std::cout << "Size of the matrix" << std::endl;
    std::cin >> n;
    HermitianMatrix matrix(n);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> dist(0, 1);
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            matrix.set_raw(i, j, dist(gen));
        }
    }

    auto start = std::chrono::high_resolution_clock::now();
    SymmetricTridiagonalMatrix stmat = matrix.tridiagonalise();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Tridiagonalisation time: " << duration.count() << " seconds." << std::endl;

    start = std::chrono::high_resolution_clock::now();
    std::vector<double> eigenvalues = stmat.eigenvalues();
    end = std::chrono::high_resolution_clock::now();
    duration = end - start;
    std::cout << "Eigenvalue time: " << duration.count() << " seconds." << std::endl;

    return 0;
}