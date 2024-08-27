#include <iostream>
#include <random>
#include <numbers>
#include <ctime>
#include <cassert>

#include "helix.h"
#include "equations.h"
#include "matrices.h"
#include "clad/Differentiator/Differentiator.h"

void Points(double points[][3], int N)
{
    // add a,b,c,d,alph, bet
    auto seed = time(nullptr);
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> uniform(-2 * std::numbers::pi_v<double>, 2 * std::numbers::pi_v<double>);
    double a = 5.2122, b = -4.79395, c = -26.40835, d = -4.207055, alph = -3.60384, bet = 1.13255;
    double t = 0;
    for (int i = 0; i < N; i++)
    {
        t += 0.1;
        double output[3];
        HelixPoint(a, b, c, d, alph, bet, t, output);
        double x = output[0], y = output[1], z = output[2];
        // double x = output[0] + (uniform(rng) / 30), y = output[1] + (uniform(rng) / 30), z = output[2] + (uniform(rng) / 30);
        points[i][0] = x;
        points[i][1] = y;
        points[i][2] = z;
        std::cout << x << " " << y << " " << z << "\n";
    }
    std::cout << "end\n";
}
double DistanceToPoint(double a, double b, double c, double d, double alph, double bet, double t, double x, double y, double z)
{
    double output[3];
    HelixPoint(a, b, c, d, alph, bet, t, output);
    return DistanceA(output, x, y, z);
}

void DistancesToAllPoints(double *points, int N, double a, double b, double c, double d, double alph, double bet, double *dist)
{
    int n = 0;
    for (int i = 0; i < N * 3; i += 3)
    {
        double x = points[i];
        double y = points[i + 1];
        double z = points[i + 2];
        double t = HelixClosestTime(a, b, c, d, alph, bet, x, y, z);
        dist[n] = DistanceToPoint(a, b, c, d, alph, bet, t, x, y, z);
        dist[n] += 0.001 * (std::abs(a) + std::abs(b) + std::abs(c) + std::abs(d) + std::abs(alph) + std::abs(bet));
        n++;
    }
}
void Jacobian(double points[][3], int N, double a, double b, double c, double d, double alph, double bet, double Jacobian[][6])
{
    /*Construct the N x 6 Jacobian.*/
    auto dist_grad = clad::gradient(DistanceToPoint, "a, b, c, d, alph, bet");

    for (int i = 0; i < N; i++)
    {
        double x = points[i][0];
        double y = points[i][1];
        double z = points[i][2];
        double t = HelixClosestTime(a, b, c, d, alph, bet, x, y, z);
        double da = 0, db = 0, dc = 0, dd = 0, dalph = 0, dbet = 0;
        dist_grad.execute(a, b, c, d, alph, bet, t, x, y, z, &da, &db, &dc, &dd, &dalph, &dbet);
        Jacobian[i][0] = da;
        Jacobian[i][1] = db;
        Jacobian[i][2] = dc;
        Jacobian[i][3] = dd;
        Jacobian[i][4] = dalph;
        Jacobian[i][5] = dbet;
    }
}

double LambdaChange(double *points, int N, double &a, double &b, double &c, double &d, double &alph, double &bet, double &lambda, double &old_square_err, double *results)
{
    double dist[N];
    DistancesToAllPoints(points, N, a + results[0], b + results[1], c + results[2], d + results[3], alph + results[4], bet + results[5], dist);
    double square_err = 0;
    double new_lambda;
    for (int i = 0; i < N; i++)
    {
        square_err += (dist[i] * dist[i]);
    }
    std::cerr << "SQUARE ERR " << square_err << std::endl;
    if (square_err >= old_square_err)
        new_lambda = lambda * 2;
    else
    {
        std::cerr << "IMPROVEMENTS!";

        a += results[0];
        b += results[1];
        c += results[2];
        d += results[3];
        alph += results[4];
        bet += results[5];
        new_lambda = lambda / 2;
        old_square_err = square_err;
    }
    std::cerr << "CHANGE IN RESULTS" << std::endl;
    for (int i = 0; i < 6; i++)
    {
        std::cerr << results[i] << " ";
    }
    std::cerr << std::endl;

    std::cerr << "RESULTS! " << std::endl;
    std::cerr << "RESULTS! " << a << " " << b << " " << c << " " << d << " " << alph << " " << bet << " " << std::endl;
    lambda = new_lambda;
    return new_lambda - lambda;
}
void LevenbergMarquardt(double points[][3], int N)
{
    double points1D[N * 3];
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            points1D[i * 3 + j] = points[i][j];
        }
    }
    auto seed = time(nullptr);
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> uniform(-2 * std::numbers::pi_v<double>, 2 * std::numbers::pi_v<double>);
    // double a = uniform(rng), b = uniform(rng), c = uniform(rng), d = uniform(rng), alph = uniform(rng), bet = uniform(rng);
    double a = 8.2122, b = -4.79395, c = -16.40835, d = -8.207055, alph = -2.60384, bet = 2.13255;
    // double a = 5.2122, b = -4.79395, c = -26, d = -4.207055, alph = -3.60384, bet = 1.13255;

    double lambda = 1;
    double lambda_change = 1;
    double old_square_err;
    double jacobian[N][6];
    double jacobian1D[N * 6];
    double tjacobian[6 * N];
    double tjj[6 * 6];
    double results[6];
    {
        double dist[N];
        DistancesToAllPoints(points1D, N, a, b, c, d, alph, bet, dist);
        old_square_err = 0;
        for (int i = 0; i < N; i++)
        {
            old_square_err += (dist[i] * dist[i]);
        }
    }
    //  --------------------------------------------------------------------------------------------------------------------
    int nr = 0;
    while (true)
    {

        if ((std::abs(lambda_change) < 0.05 && lambda_change < 0) || nr > 2000)
            break;

        Jacobian(points, N, a, b, c, d, alph, bet, jacobian);
        std::cerr << "Jacobian: " << std::endl;
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < 6; j++)
            {
                jacobian1D[i * 6 + j] = jacobian[i][j];
                std::cerr << jacobian[i][j] << " ";
            }
            std::cerr << std::endl;
        }
        Transpose(jacobian1D, N, 6, tjacobian);

        MatrixMultiply(tjacobian, jacobian1D, 6, N, 6, tjj);
        double diag[6 * 6];
        DiagOfSquareM(tjj, 6, diag);

        double identity[6 * 6];
        ScalarMultiply(diag, 6, 6, lambda, identity);
        double left_side[6 * 6];
        AddMatrices(tjj, identity, 6, 6, left_side);
        double dist[N];
        DistancesToAllPoints(points1D, N, a, b, c, d, alph, bet, dist);
        double right_side[6 * 1];
        std::cerr << "dist " << std::endl;

        for (int j = 0; j < N; j++)
        {

            std::cerr << dist[j] << "  ";
        }
        std::cerr << std::endl;
        MatrixMultiply(tjacobian, dist, 6, N, 1, right_side);
        ScalarMultiply(right_side, 1, 6, -1, right_side);

        // left side is 6x6, right side is 6x1, so h is 6x1.
        double left_side3D[6][6];
        std::cerr << "lambda " << lambda << std::endl;

        std::cerr << "left_side3D " << std::endl;
        for (int i = 0; i < 6; i++)
        {
            for (int j = 0; j < 6; j++)
            {
                left_side3D[i][j] = left_side[i * 6 + j];
                std::cerr << left_side3D[i][j] << "  ";
            }
            std::cerr << std::endl;
        }
        double forward_elim[6][6];
        double unchanged_rs[6];
        CopyMatrix(right_side, 6, unchanged_rs);
        ForwardElim(left_side3D, right_side, forward_elim);
        BackSub(forward_elim, right_side, results);
        CheckSolution(left_side3D, unchanged_rs, results);
        lambda_change = LambdaChange(points1D, N, a, b, c, d, alph, bet, lambda, old_square_err, results);
        nr++;
    }
    double t = -N;
    for (int i = 0; i < 20 * N; i++)
    {
        t += 0.1;
        double output[3];
        HelixPoint(a, b, c, d, alph, bet, t, output);
        double x = output[0], y = output[1], z = output[2];

        std::cout << x << " " << y << " " << z << "\n";
    }
    std::cout << "end\n";
}
