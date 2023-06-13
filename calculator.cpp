#include <iostream>
#include <vector>

class Vector3D {
private:
    std::vector<double> coordinates;

public:
    Vector3D(double x, double y, double z) : coordinates({x, y, z}) {
    }

    double& operator[](int index) {
        return coordinates[index];
    }

    const double& operator[](int index) const {
        return coordinates[index];
    }

    Vector3D operator+(const Vector3D& other) const {
        double x = coordinates[0] + other.coordinates[0];
        double y = coordinates[1] + other.coordinates[1];
        double z = coordinates[2] + other.coordinates[2];
        return Vector3D(x, y, z);
    }

    Vector3D operator-(const Vector3D& other) const {
        double x = coordinates[0] - other.coordinates[0];
        double y = coordinates[1] - other.coordinates[1];
        double z = coordinates[2] - other.coordinates[2];
        return Vector3D(x, y, z);
    }

    Vector3D operator*(double scalar) const {
        double x = coordinates[0] * scalar;
        double y = coordinates[1] * scalar;
        double z = coordinates[2] * scalar;
        return Vector3D(x, y, z);
    }

    Vector3D operator/(double scalar) const {
        double x = coordinates[0] / scalar;
        double y = coordinates[1] / scalar;
        double z = coordinates[2] / scalar;
        return Vector3D(x, y, z);
    }

    friend std::ostream& operator<<(std::ostream& output, const Vector3D& vector) {
        output << "(" << vector[0] << ", " << vector[1] << ", " << vector[2] << ")";
        return output;
    }
};

class Matrix {
private:
    std::vector<std::vector<double>> elements;

public:
    Matrix(const std::vector<std::vector<double>>& values) : elements(values) {
    }

    std::vector<double>& operator()(int row) {
        return elements[row];
    }

    const std::vector<double>& operator()(int row) const {
        return elements[row];
    }

    int numRows() const {
        return elements.size();
    }

    int numCols() const {
        if (numRows() > 0) {
            return elements[0].size();
        }
        return 0;
    }

    Matrix operator+(const Matrix& other) const {
        if (numRows() != other.numRows() || numCols() != other.numCols()) {
            throw std::runtime_error("Matrix dimensions do not match.");
        }

        std::vector<std::vector<double>> result(numRows(), std::vector<double>(numCols(), 0.0));

        for (int i = 0; i < numRows(); i++) {
            for (int j = 0; j < numCols(); j++) {
                result[i][j] = elements[i][j] + other.elements[i][j];
            }
        }

        return Matrix(result);
    }

    Matrix operator-(const Matrix& other) const {
        if (numRows() != other.numRows() || numCols() != other.numCols()) {
            throw std::runtime_error("Matrix dimensions do not match.");
        }

        std::vector<std::vector<double>> result(numRows(), std::vector<double>(numCols(), 0.0));

        for (int i = 0; i < numRows(); i++) {
            for (int j = 0; j < numCols(); j++) {
                result[i][j] = elements[i][j] - other.elements[i][j];
            }
        }

        return Matrix(result);
    }

    Matrix operator*(double scalar) const {
        std::vector<std::vector<double>> result(numRows(), std::vector<double>(numCols(), 0.0));

        for (int i = 0; i < numRows(); i++) {
            for (int j = 0; j < numCols(); j++) {
                result[i][j] = elements[i][j] * scalar;
            }
        }

        return Matrix(result);
    }

    Matrix operator/(double scalar) const {
        if (scalar == 0.0) {
            throw std::runtime_error("Division by zero.");
        }

        std::vector<std::vector<double>> result(numRows(), std::vector<double>(numCols(), 0.0));

        for (int i = 0; i < numRows(); i++) {
            for (int j = 0; j < numCols(); j++) {
                result[i][j] = elements[i][j] / scalar;
            }
        }

        return Matrix(result);
    }

    friend std::ostream& operator<<(std::ostream& output, const Matrix& matrix) {
        for (int i = 0; i < matrix.numRows(); i++) {
            for (int j = 0; j < matrix.numCols(); j++) {
                output << matrix.elements[i][j] << " ";
            }
            output << std::endl;
        }
        return output;
    }
};

int main() {
    Vector3D v1(1.0, 2.0, 3.0);
    Vector3D v2(4.0, 5.0, 6.0);

    std::cout << "v1: " << v1 << std::endl;
    std::cout << "v2: " << v2 << std::endl;

    Vector3D sum = v1 + v2;
    std::cout << "v1 + v2: " << sum << std::endl;

    Vector3D difference = v1 - v2;
    std::cout << "v1 - v2: " << difference << std::endl;

    Vector3D scaled = v1 * 2.0;
    std::cout << "v1 * 2.0: " << scaled << std::endl;

    Vector3D divided = v2 / 3.0;
    std::cout << "v2 / 3.0: " << divided << std::endl;

    std::cout << std::endl;

    Matrix m1({{1.0, 2.0}, {3.0, 4.0}});
    Matrix m2({{5.0, 6.0}, {7.0, 8.0}});

    std::cout << "m1:" << std::endl;
    std::cout << m1 << std::endl;

    std::cout << "m2:" << std::endl;
    std::cout << m2 << std::endl;

    Matrix sumMatrix = m1 + m2;
    std::cout << "m1 + m2:" << std::endl;
    std::cout << sumMatrix << std::endl;

    Matrix diffMatrix = m1 - m2;
    std::cout << "m1 - m2:" << std::endl;
    std::cout << diffMatrix << std::endl;

    Matrix scaledMatrix = m1 * 2.0;
    std::cout << "m1 * 2.0:" << std::endl;
    std::cout << scaledMatrix << std::endl;

    Matrix dividedMatrix = m2 / 2.0;
    std::cout << "m2 / 2.0:" << std::endl;
    std::cout << dividedMatrix << std::endl;

    return 0;
}