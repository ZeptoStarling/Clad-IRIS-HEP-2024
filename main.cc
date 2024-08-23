#include <iostream>
#include <random>
#include <numbers>
#include <ctime>

#include "helix.h"
#include "equations.h"
#include "matrices.h"

void Points(double points[][3], int N)
{
	auto seed = time(nullptr);
	std::mt19937_64 rng(seed);
	std::uniform_real_distribution<double> uniform(-2 * std::numbers::pi_v<double>, 2 * std::numbers::pi_v<double>);
	double a = uniform(rng), b = uniform(rng), c = uniform(rng), d = uniform(rng), alph = uniform(rng), bet = uniform(rng);
	double t = 0;
	for (int i = 0; i < N; i++)
	{
		t += 0.1;
		double output[3];
		HelixPoint(a, b, c, d, alph, bet, t, output);
		double x = output[0] + (uniform(rng) / 20), y = output[1] + (uniform(rng) / 20), z = output[2] + (uniform(rng) / 20);
		points[i][0] = x;
		points[i][1] = y;
		points[i][2] = z;
	}
}
double DistanceToPoint(double a, double b, double c, double d, double alph, double bet, double t, double x, double y, double z)
{
	double output[3];
	HelixPoint(a, b, c, d, alph, bet, t, output);
	return DistanceSquareA(output, x, y, z);
}
void Jacobian(double points[][3], int N, double a, double b, double c, double d, double alph, double bet, double Jacobian[][6])
{
	/*Construct the Jacobian.*/
	auto dist_grad = clad::gradient(DistanceToPoint, "a, b, c, d, alph, bet");

	for (int i = 0; i < N; i++)
	{
		double x = points[i][0];
		double y = points[i][1];
		double z = points[i][2];
		double t = HelixClosestTime(a, b, c, d, alph, bet, x, y, z);
		double distFound = DistanceToPoint(a, b, c, d, alph, bet, t, x, y, z);
		double da = 0, db = 0, dc = 0, dd = 0, dalph = 0, dbet = 0;
		dist_grad.execute(a, b, c, d, alph, bet, t, x, y, z, &da, &db, &dc, &dd, &dalph, &bet);
		Jacobian[i][0] = da;
		Jacobian[i][1] = db;
		Jacobian[i][2] = dc;
		Jacobian[i][3] = dd;
		Jacobian[i][4] = dalph;
		Jacobian[i][5] = dbet;
	}
}
void LevenbergMarquardt(double points[][3], int N)
{
	// TO DO. Fix issues with matrices.h, check if output is correct
	// TO DO. Not finished. No solution for the actual equation, no updating of lambda.
	auto seed = time(nullptr);
	std::mt19937_64 rng(seed);
	std::uniform_real_distribution<double> uniform(-2 * std::numbers::pi_v<double>, 2 * std::numbers::pi_v<double>);
	double a = uniform(rng), b = uniform(rng), c = uniform(rng), d = uniform(rng), alph = uniform(rng), bet = uniform(rng);
	double jacobian[N][6];
	Jacobian(points, N, a, b, c, d, alph, bet, jacobian);
	double tjacobian[6][200];
	Transpose(jacobian, tjacobian);

	double lambda = 0.001;
	double tjj[6][6];
	for (int i = 0; i < 6; i++)
	{
		for (int j = 0; j < 6; j++)
		{
			tjj[i][j] = 0.0;
			for (int k = 0; k < 6; k++)
			{
				tjj[i][j] += tjacobian[i][k] * jacobian[k][j];
			}
		}
	}
	// MatrixMultiply(tjacobian, jacobian, tjj);
	double diag[6][6];
	Diag(tjj, diag);
	double left_side[6][6];
	ScalarMultiply(diag, lambda, left_side);
}
int main()
{
	int N = 200;
	double points[N][3];
	Points(points, N);
	LevenbergMarquardt(points, N);
}
