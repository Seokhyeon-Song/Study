#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <algorithm>
#include <array>
#include <complex>
#include <cstdlib>
#include <string>

constexpr auto POINTS = 200; // Starts from 0, ends at POINTS
constexpr auto SEARCH_POINTS = 10; // -SEARCH_POINTS to SEARCH_POINTS, for each of the two parameters
constexpr auto SIGMA_PRECISION = 4;
constexpr auto CONVOLUTION_POINTS = 100;
constexpr auto ITERS = 5;
constexpr auto PLOT_POINTS = 100;

using namespace std::complex_literals;
using comp = std::complex<double>;

std::array<double, POINTS + 1> pValues, resonanceData, narrowWidthData;
std::array<double, CONVOLUTION_POINTS + 1> convolutionValues;

double massOS, widthOS, stDev, gaussianNorm;

double LinearCoeff(const std::array<double, POINTS + 1>& X, const std::array<double, POINTS + 1>& Y) {
	double num = 0, den = 0;
	for (int i = 0; i <= POINTS; i++) {
		num += X[i] * Y[i];
		den += X[i] * X[i];
	}
	return num / den;
}

inline bool FloatEqual(double a, double b) {
	return fabs(a - b) < FLT_EPSILON * std::fmax(a,b);
}

inline double NWA(double p, double a, double b) {
	return 1 / (pow(p * p - a * a, 2) - b * b);
}

comp F0(double m, double p) {
	if (FloatEqual(p, 0)) {
		return -2.0 * log(comp(m));
	}
	else if (FloatEqual(p, 2 * m)) {
		return 2.0 - 2.0 * log(comp(m));
	}
	else {
		return 2.0 - 2.0 * log(comp(m)) + (2.0 * sqrt(comp(4.0 * m * m - p * p)) / p) * atan(-p / sqrt(comp(4.0 * m * m - p * p)));
	}
}

comp FPrime0(double m, double p) {
	if (FloatEqual(p, 2 * m)) {
		return comp(-1 / (4 * m * m));
	}
	else {
		return -1.0 / comp(p * p) - (4.0 * m * m / (sqrt(comp(4.0 * m * m - p * p)) * p * p * p)) * atan(-p / sqrt(comp(4.0 * m * m - p * p)));
	}
}

comp PropagatorScalar(double p, double M, double m, double c, double gam) {
	comp den = p * p - M * M + 1.0i * M * gam + c * (F0(m, p) - F0(m, M) - (p * p - M * M) * real(FPrime0(m, M)));
	return 1.0 / den;
}

double ResonanceConvolution(double p, double M, double m, double c, double gam) {
	double res = 0.0;
	for (int i = -CONVOLUTION_POINTS; i <= CONVOLUTION_POINTS; i++) {
		res += convolutionValues[abs(i)] * pow(abs(PropagatorScalar(p + i * stDev * SIGMA_PRECISION / CONVOLUTION_POINTS, M, m, c, gam)), 2);
	}
	return res * gaussianNorm;
}

double NWAConvolution(double p, double a, double b) {
	double res = 0.0;
	for (int i = -CONVOLUTION_POINTS; i <= CONVOLUTION_POINTS; i++) {
		res += convolutionValues[abs(i)] * NWA(p + i * stDev * SIGMA_PRECISION / CONVOLUTION_POINTS, a, b);
	}
	return res * gaussianNorm;
}

double diffSqSum(const std::array<double, POINTS + 1>& dat, double a, double b) {
	for (int i = 0; i <= POINTS; i++) {
		narrowWidthData[i] = NWAConvolution(pValues[i], a, b);
	}

	double cCurrent = LinearCoeff(narrowWidthData, resonanceData);
	double res = 0.0;
	for (int i = 0; i <= POINTS; i++) {
		res += pow(resonanceData[i] - cCurrent * narrowWidthData[i], 2);
	}
	return res;
}

double diffSectionRatio(double massChi, double coup) {
	for (int i = 0; i <= POINTS; i++) {
		resonanceData[i] = ResonanceConvolution(pValues[i], massOS, massChi, coup, widthOS);
	}

	double aCentre = massOS, bCentre = massOS * widthOS;
	double aStep = 2 * widthOS / SEARCH_POINTS, bStep = massOS * widthOS * SEARCH_POINTS / ((SEARCH_POINTS + 1) * SEARCH_POINTS);
	double aCurrent, bCurrent;
	int minElem;
	std::array<double, (2 * SEARCH_POINTS + 1)* (2 * SEARCH_POINTS + 1)> cost;

	for (int k = 0; k < ITERS; k++) {
		for (int i = -SEARCH_POINTS; i <= SEARCH_POINTS; i++) {
			aCurrent = aCentre + aStep * i;
			for (int j = -SEARCH_POINTS; j <= SEARCH_POINTS; j++) {
				bCurrent = bCentre + bStep * j;
				cost[(i + SEARCH_POINTS) * (2 * SEARCH_POINTS + 1) + j + SEARCH_POINTS] = diffSqSum(resonanceData, aCurrent, bCurrent);
			}
		}
		minElem = std::min_element(cost.begin(), cost.end()) - cost.begin();
		aCentre = aCentre + aStep * (minElem / (2 * SEARCH_POINTS + 1) - SEARCH_POINTS);
		bCentre = bCentre + bStep * (minElem % (2 * SEARCH_POINTS + 1) - SEARCH_POINTS);
		aStep /= SEARCH_POINTS;
		bStep /= SEARCH_POINTS;
	}

	for (int i = 0; i <= POINTS; i++) {
		narrowWidthData[i] = NWAConvolution(pValues[i], aCentre, bCentre);
	}
	double cCurrent = LinearCoeff(narrowWidthData, resonanceData);

	double resonanceTotal = 0.0, diffTotal = 0.0;
	for (int i = 0; i <= POINTS; i++) {
		resonanceTotal += resonanceData[i];
		diffTotal += abs(resonanceData[i] - cCurrent * narrowWidthData[i]);
	}
	return diffTotal / resonanceTotal;
}

int main(void) {
	massOS = 400, widthOS = 40, stDev = 20;
	gaussianNorm = SIGMA_PRECISION / (CONVOLUTION_POINTS * sqrt(2 * M_PI));
	for (int i = 0; i <= CONVOLUTION_POINTS; i++) {
		convolutionValues[i] = exp(-pow((double)i * SIGMA_PRECISION / CONVOLUTION_POINTS, 2) / 2);
	}
	for (int i = 0; i <= POINTS; i++) {
		pValues[i] = massOS + 8 * widthOS * ((double)i / POINTS - 0.5);
	}
	
	std::string filename = "fitresult.csv";
	std::ofstream writeFile;
	writeFile.open(filename);
	writeFile << "m, c, val," << std::endl;
	int totalPoints = PLOT_POINTS * PLOT_POINTS;
	int percent = 0;
	for (int i = 0; i < PLOT_POINTS; i++) {
		double massChi = massOS / 2 + widthOS * (i - PLOT_POINTS / 2) / PLOT_POINTS;
		for (int j = 1; j <= PLOT_POINTS; j++) {
			double coup = massOS * massOS * j / (4 * PLOT_POINTS * M_PI);
			writeFile << massChi << ", " << coup << ", " << diffSectionRatio(massChi, coup) << "," << std::endl;
			if (percent < (PLOT_POINTS * i + j) * 100 / totalPoints) {
				percent = (PLOT_POINTS * i + j) * 100 / totalPoints;
				std::cout << percent << "% Done" << std::endl;
			}
		}
	}
	writeFile.close();
	return 0;
}