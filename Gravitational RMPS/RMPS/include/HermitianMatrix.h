#ifndef HERMITIAN_MATRIX_H
#define HERMITIAN_MATRIX_H

#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <limits>

using complex = std::complex<double>;

complex dotProduct(const std::vector<complex>& a, const std::vector<complex>& b) {
    complex res;
    // Caution: No size check for efficiency.
    for (size_t i = 0; i < a.size(); i++) {
        res += conj(a[i]) * b[i];
    }
    return res;
}

class SymmetricTridiagonalMatrix {
public:
    std::vector<double> diag;  // Main diagonal (size n)
    std::vector<double> subdiag; // Subdiagonal (size n-1)
    size_t size;

    // Constructor
    SymmetricTridiagonalMatrix(size_t n) : size(n), diag(n), subdiag(n - 1) {}

    // Get full matrix element (i, j)
    double getElement(size_t i, size_t j) const {
        if (i >= diag.size() || j >= diag.size()) throw std::out_of_range("Index out of bounds");
        if (i == j) return diag[i]; // Diagonal elements
        if (i + 1 == j || j + 1 == i) return subdiag[std::min(i, j)]; // Sub/Superdiagonal elements
        return 0.0; // Zero elsewhere
    }

    // Print the matrix (for debugging)
    void print() const {
        size_t n = diag.size();
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                std::cout << getElement(i, j) << " ";
            }
            std::cout << "\n";
        }
    }

    void iterateQR();

    std::vector<double> eigenvalues();
};

class HermitianMatrix {
public:
    std::vector<double> data;
    size_t size;

    // Constructor: initializes a Hermitian matrix of a given size
    explicit HermitianMatrix(size_t n) : size(n), data(n * n) {}

    // Accessor function to get the size of the matrix
    size_t getSize() const { return size; }

    // Function to set an element
    void set(size_t i, size_t j, const complex& value) {
        if (i >= size || j >= size) {
            throw std::out_of_range("Index out of range");
        }
        if (i == j) {
            data[i * size + i] = real(value);
        } else if (i > j) {
            data[i * size + j] = real(value);
            data[j * size + i] = imag(value);
        } else {
            data[j * size + i] = real(value);
            data[i * size + j] = -imag(value);
        }
    }

    void set_raw(size_t i, size_t j, const double value) {
        data[i * size + j] = value;
    }

    // Function to get an element from the matrix
    complex get(size_t i, size_t j) const {
        if (i >= size || j >= size) {
            throw std::out_of_range("Index out of range");
        }
        if (i == j) {
            return complex(data[i * size + i], 0.0);
        } else if (i > j) {
            return complex(data[i * size + j], data[j * size + i]);
        } else if (i < j) {
            return complex(data[j * size + i], -data[i * size + j]);
        }
    }

    double get_raw(size_t i, size_t j) const {
        return data[i * size + j];
    }

    // Function to print the Hermitian matrix
    void print() const {
        for (size_t i = 0; i < size; ++i) {
            for (size_t j = 0; j < size; ++j) {
                std::cout << get(i, j) << " ";
            }
            std::cout << std::endl;
        }
    }

    HermitianMatrix operator*(const HermitianMatrix& B) const {
        HermitianMatrix result(size);
        for (size_t i = 0; i < size; ++i) {
            for (size_t j = 0; j <= i; ++j) {
                complex sum = 0;
                for (size_t k = 0; k < size; ++k) {
                    sum += get(i, k) * B.get(k, j);
                }
                result.set(i, j, sum);
            }
        }
        return result;
    }

    SymmetricTridiagonalMatrix tridiagonalise() {
        SymmetricTridiagonalMatrix res(1);
        res.size = size;
        res.diag[0] = data[0];
        for (size_t k = 0; k < size - 2; k++) {
            // Computing Householder vector x
            std::vector<complex> x;
            double xNorm = 0.0;
            for (size_t j = k + 1; j < size; j++) {
                x.push_back(complex(data[j*size+k], data[k*size+j]));
                xNorm += data[j*size+k] * data[j*size+k] + data[k*size+j] * data[k*size+j];
            }
            double sigma = std::sqrt(xNorm);
            xNorm -= norm(x[0]);
            x[0] -= sigma;
            xNorm += norm(x[0]);
            double xInvMag = 1 / std::sqrt(xNorm);
            for (size_t i = 0; i < x.size(); i++) {
                x[i] *= xInvMag;
            }

            // Applying Householder transformation
            std::vector<complex> y; // (y0, y)=Ax
            complex y0;
            for (size_t j = 0; j < x.size(); j++) y0 += complex(get_raw(j + k + 1, k), -get_raw(k, j + k + 1)) * x[j];
            for (size_t i = k + 1; i < size; i++) {
                complex t;
                for (size_t j = 0; j < x.size(); j++) {
                    t += get(i, j + k + 1) * x[j];
                }
                y.push_back(t);
            }
            complex xdoty = dotProduct(x, y);
            for (size_t i = k + 1; i < size; i++) {
                for (size_t j = i; j < size; j++) {
                    set(i, j, get(i, j) - 2.0 * x[i-k-1] * conj(y[j-k-1]) - 2.0 * y[i-k-1] * conj(x[j-k-1]) + 4.0 * x[i-k-1] * conj(x[j-k-1]) * xdoty);
                }
            }
            res.subdiag.push_back(abs(get(k, k+1) - 2.0 * y0 * conj(x[0])));
            res.diag.push_back(data[(k + 1) * (size + 1)]);
        }
        res.subdiag.push_back(abs(get(size - 2, size - 1)));
        res.diag.push_back(data[size * size - 1]);
        return res;
    }
};

void SymmetricTridiagonalMatrix::iterateQR() {
    // Caution: Assuming that the matrix is tridiagonal, but no check for efficiency.
    // Choosing the shift
    double mu;
    if (diag[size-1] > diag[size-2]) {
        mu = diag[size-1] + diag[size-2] + std::sqrt((diag[size-1] - diag[size-2]) * (diag[size-1] - diag[size-2]) + 4 * subdiag[size-2] * subdiag[size-2]);
        mu /= 2;
    } else {
        mu = diag[size-1] + diag[size-2] - std::sqrt((diag[size-1] - diag[size-2]) * (diag[size-1] - diag[size-2]) + 4 * subdiag[size-2] * subdiag[size-2]);
        mu /= 2;
    }

    // First Givens rotation
    double r = std::sqrt((diag[0] - mu) * (diag[0] - mu) + subdiag[0] * subdiag[0]);
    double bulge;
    if (r < std::numeric_limits<double>::epsilon()) std::cout<< "Divide by 0 error during rotation" << std::endl;
    double c = (diag[0] - mu) / r;
    double s = subdiag[0] / r;
    diag[0] += (1 - c * c) * (diag[0] + diag[1] - 2 * mu);
    subdiag[0] = c * (diag[1] - mu) * s + (c * c - 1) * subdiag[0];
    diag[1] = c * c * (diag[1] - mu) + (c * c - 1) * r * c + mu;
    bulge = s * subdiag[1];
    subdiag[1] *= c;

    // Bulge chasing
    for (size_t k = 0; k < size - 2; k++) {
        r = std::sqrt(subdiag[k] * subdiag[k] + bulge * bulge);
        if (r < std::numeric_limits<double>::epsilon()) std::cout<< "Divide by 0 error during chasing" << std::endl;
        double c = subdiag[k] / r;
        double s = bulge / r;
        subdiag[k] = r;
        if (k < size - 3) {
            bulge = s * subdiag[k+2];
            subdiag[k+2] *= c;
        }
        double a = diag[k+1, k+1], d = diag[k+2];
        double b = subdiag[k+1];
        diag[k+1] = c * c * a + s * s * d + 2 * b * c * s;
        diag[k+2] = s * s * a + c * c * d - 2 * b * c * s;
        subdiag[k+1] = c * s * (d - a) + (2 * c * c - 1) * b;
    }
}

std::vector<double> SymmetricTridiagonalMatrix::eigenvalues() {
    std::vector<double> res;
    for (; size > 2; size--) {
        do { iterateQR(); } while (abs(subdiag[size-2]) > std::numeric_limits<float>::epsilon());
        res.push_back(diag[size-1]);
        diag.pop_back();
        subdiag.pop_back();
    }
    res.push_back(0.5 * (diag[1] + diag[0] + std::sqrt((diag[1] - diag[0]) * (diag[1] - diag[0]) + 4 * subdiag[0] * subdiag[0])));
    res.push_back(0.5 * (diag[1] + diag[0] - std::sqrt((diag[1] - diag[0]) * (diag[1] - diag[0]) + 4 * subdiag[0] * subdiag[0])));
    return res;
}

#endif // HERMITIAN_MATRIX_H