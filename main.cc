#include <iostream>
#include <cmath>
#include "clad/Differentiator/Differentiator.h"

#include "fitter.h"

int main()
{
	int nr_of_points = 200;
	double points[nr_of_points * 3];
	Points(points, nr_of_points);
	// LevenbergMarquardt(points, nr_of_points);
	GradientDescent(points, nr_of_points);
	// for (int i = 0; i < nr_of_points; i++)
	// {
	// 	std::cout << points[i * 3 + 0] << " " << points[i * 3 + 1] << " " << points[i * 3 + 2] << "\n";
	// }
	// std::cout << "end\n";
	double a = 5.2122, b = -4.79395, c = -26.40835, d = -4.207055, alph = -3.60384, bet = 1.13255;
	// GridSearchParams(points, nr_of_points, a, b, c, d, alph, bet);
	// double dist[nr_of_points];
	// DistancesToAllPoints(points, nr_of_points, a, 0.1, c, d, alph, bet, dist);
	// PrintMatrix("dist ", dist, 200, 1);
}
