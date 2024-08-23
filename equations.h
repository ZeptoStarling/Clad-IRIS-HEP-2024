#pragma once

#include <cmath>
#include <cassert>

#include "clad/Differentiator/Differentiator.h"

double EvaluateSinPlusLin(double A, double B, double C, double D, double x) {
	return A * std::sin(x + B) + C * x + D;
}

void PrintSinPlusLin(double A, double B, double C, double D) {
	std::cerr << A << " * sin(x + " << B << ") + " << C << " * x + " << D << " = 0\n";
}

// A sin (x + B) + C x + D = 0
double SolveSinPlusLin(double A, double B, double C, double D, double mi = -1000000, double ma = 1000000) {
	for (int i = 0; i < 100; i++) {
		double mid = (mi + ma) / 2;
		double vmi = EvaluateSinPlusLin(A, B, C, D, mi);
		double vmid = EvaluateSinPlusLin(A, B, C, D, mid);
		double vma = EvaluateSinPlusLin(A, B, C, D, ma);

		if (vmi < 0 and 0 < vmid) {
			ma = mid;
		}
		else if (vmid < 0 and 0 < vma) {
			mi = mid;
		}
		else if (vmid < 0 and 0 < vmi) {
			ma = mid;
		}
		else if (vma < 0 and 0 < vmid) {
			mi = mid;
		}
		else {
			break;
			mi = mid;
		}
	}

	double x = (mi + ma) / 2;
	return x;
	auto derivative = clad::differentiate(EvaluateSinPlusLin, "x");

	for (int i = 0; i < 100; i++) {
		x = x - EvaluateSinPlusLin(A, B, C, D, x) / derivative.execute(A, B, C, D, x);
	}

	// assert(std::abs(EvaluateSinPlusLin(A, B, C, D, x)) < 1e-10);
	return x;
}

double NextValPiK(double offs, double x) {
	constexpr auto MY_PI = std::numbers::pi_v<double>;

	if (x < 0) {
		double v = -NextValPiK(-offs, -x)+ 2 * MY_PI;
		return v > x ? v : v + 2 * MY_PI;
	}

	double kie = std::floor(x / 2 / MY_PI);

	for (int i = -2; i <= 2; i++) {
		double v = (kie + i) * 2 * MY_PI + offs;

		if (v > x) {
			return v;
		}
	}

	return 1000000000;
}

double PrevValPiK(double offs, double x) {
	constexpr auto MY_PI = std::numbers::pi_v<double>;
	return NextValPiK(offs, x) - 2 * MY_PI;
}

// A cos(x + B) + C = 0
double NextSinPlusInflection(double A, double B, double C, double x) {
	// cos(x + B) = -C / A
	double inv = std::acos(-C / A);
	return std::min(NextValPiK(inv - B, x), NextValPiK(-inv - B, x));
}

// A cos(x + B) + C = 0
double PrevSinPlusInflection(double A, double B, double C, double x) {
	// cos(x + B) = -C / A
	double inv = std::acos(-C / A);
	return std::max(PrevValPiK(inv - B, x), PrevValPiK(-inv - B, x));
}
