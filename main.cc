#include <iostream>
#include <cmath>
#include "clad/Differentiator/Differentiator.h"

#include "fitter.h"

int main()
{
	int N = 10;
	double points[N][3];
	Points(points, N);
	for (int i = 0; i < N; i++)
	{
		std::cerr << points[i][0] << " " << points[i][1] << " " << points[i][2] << std::endl;
	}
	LevenbergMarquardt(points, N);
}
