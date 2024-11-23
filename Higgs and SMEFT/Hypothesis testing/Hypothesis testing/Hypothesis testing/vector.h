#pragma once
class Vector {
public:
    static const int dim = 3;

    // Vector Public Methods
    Vector() : x(0), y(0), z(0) {}
    Vector(double x, double y, double z) : x(x), y(y), z(z) {}

    double operator[](int i) const {
        if (i == 0) return x;
        if (i == 1) return y;
        return z;
    }
    double& operator[](int i) {
        if (i == 0) return x;
        if (i == 1) return y;
        return z;
    }

    Vector operator+(Vector c) const {
        return Vector(x + c.x, y + c.y, z + c.z);
    }
    Vector operator-(Vector c) const {
        return Vector(x - c.x, y - c.y, z - c.z);
    }

    // Vector Public Members
    double x, y, z;
};

class Matrix {
public:
    Matrix() : matrix{ {1,0,0},{0,1,0},{0,0,1} }, invmat{ {1,0,0},{0,1,0},{0,0,1} } {}

    Matrix(const double m[3][3]) {
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                matrix[i][j] = m[i][j];
        invmat[0][0] = m[1][1] * m[2][2] - m[2][1] * m[1][2];
        invmat[0][1] = m[2][1] * m[0][2] - m[0][1] * m[2][2];
        invmat[0][2] = m[0][1] * m[1][2] - m[1][1] * m[0][2];
        double invdet = 1.0 /
            (invmat[0][0] * m[0][0] + invmat[0][1] * m[1][0] + invmat[0][2] * m[2][0]);
        invmat[0][0] *= invdet;
        invmat[0][1] *= invdet;
        invmat[0][2] *= invdet;
        invmat[1][0] = (m[1][2] * m[2][0] - m[2][2] * m[1][0]) * invdet;
        invmat[1][1] = (m[2][2] * m[0][0] - m[0][2] * m[2][0]) * invdet;
        invmat[1][2] = (m[0][2] * m[1][0] - m[1][2] * m[0][0]) * invdet;
        invmat[2][0] = (m[1][0] * m[2][1] - m[2][0] * m[1][1]) * invdet;
        invmat[2][1] = (m[2][0] * m[0][1] - m[0][0] * m[2][1]) * invdet;
        invmat[2][2] = (m[0][0] * m[1][1] - m[1][0] * m[0][1]) * invdet;
    }

    Matrix(const double m[3][3], const double minv[3][3]) {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                matrix[i][j] = m[i][j];
                invmat[i][j] = minv[i][j];
            }
        }
    }

    // Transformation Public Members
    double matrix[3][3];
    double invmat[3][3];
};

Vector Transform(Vector v, Matrix t) {
    Vector res;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            res[i] += t.matrix[i][j] * v[j];
    return res;
}

Vector operator*(double c, const Vector v) {
    return Vector(c * v.x, c * v.y, c * v.z);
}