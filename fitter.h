#include <iostream>
#include <random>
#include <numbers>
#include <ctime>
#include <cassert>

#include "helix.h"
#include "equations.h"
#include "matrices.h"
#include "clad/Differentiator/Differentiator.h"
//---------------------
#include <fstream>
double SquareErr(double *points, int nr_of_points, double a, double b, double c, double d, double alph, double bet);
void Points(double *points, int nr_of_points)
{
    // add a,b,c,d,alph, bet
    auto seed = time(nullptr);
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> uniform(-2 * std::numbers::pi_v<double>, 2 * std::numbers::pi_v<double>);
    // double a = 5.2122, b = 2, c = 0.40835, d = 0.707055, alph = -3.60384, bet = 1.13255;
    double a = 5.2122, b = 2, c = 10.835, d = 17.07055, alph = -3.60384, bet = 1.13255; // with L-M and double a = 20.5, b = 0.001, c = 40.0, d = 15.207055, alph = -3.60384, bet = 1.13255; i get -d. why is the graph still good?
    double t = 0;
    for (int i = 0; i < nr_of_points; i++)
    {
        t += 0.1;
        double output[3];
        HelixPoint(a, b, c, d, alph, bet, t, output);
        double x = output[0], y = output[1], z = output[2];
        // double x = output[0] + (uniform(rng) / 30), y = output[1] + (uniform(rng) / 30), z = output[2] + (uniform(rng) / 30);
        points[i * 3] = x;
        points[i * 3 + 1] = y;
        points[i * 3 + 2] = z;
        // std::cout << x << " " << y << " " << z << "\n";
    }
    // std::cout << "end\n";
    std::cerr << "Target square err:\n";
    double square_err = SquareErr(points, nr_of_points, a, b, c, d, alph, bet);
    std::cerr << square_err << std::endl;
}
double DistanceToPoint(double a, double b, double c, double d, double alph, double bet, double x, double y, double z)
{
    double t = HelixClosestTime(a, b, c, d, alph, bet, x, y, z);
    double output[3];
    HelixPoint(a, b, c, d, alph, bet, t, output);
    double dist = DistanceA(output, x, y, z);
    // dist += 0.001 * ((a * a) + (b * b) + (c * c) + (d * d) + (alph * alph) + (bet * bet)) + (1 / b / b);
    dist += 0.001 * ((a * a) + (b * b) + (c * c) + (d * d) + (alph * alph) + (bet * bet));

    return dist;
}

void DistancesToAllPoints(double *points, int nr_of_points, double a, double b, double c, double d, double alph, double bet, double *dist)
{
    int n = 0;
    for (int i = 0; i < nr_of_points; i++)
    {
        double x = points[i * 3];
        double y = points[i * 3 + 1];
        double z = points[i * 3 + 2];
        dist[n] = DistanceToPoint(a, b, c, d, alph, bet, x, y, z);
        n++;
    }
}
void Jacobian(double *points, int nr_of_points, double a, double b, double c, double d, double alph, double bet, double *Jacobian)
{
    /*Construct the nr_of_points x 6 Jacobian.*/
    // auto dist_grad = clad::gradient(DistanceToPoint, "a, c, d, alph, bet");
    // dist_grad.dump();
    // dist_grad.dump();
    // double t = HelixClosestTime(a, b, c, d, alph, bet, points[0], points[1], points[2]);
    // auto dist_time = clad::gradient(NextSinPlusInflection, "A, B");
    // dist_time.dump();
    auto dist_a = clad::differentiate(DistanceToPoint, "a");
    auto dist_b = clad::differentiate(DistanceToPoint, "b");
    auto dist_c = clad::differentiate(DistanceToPoint, "c");
    auto dist_d = clad::differentiate(DistanceToPoint, "d");
    auto dist_alph = clad::differentiate(DistanceToPoint, "alph");
    auto dist_bet = clad::differentiate(DistanceToPoint, "bet");

    for (int i = 0; i < nr_of_points; i++)
    {
        double x = points[i * 3];
        double y = points[i * 3 + 1];
        double z = points[i * 3 + 2];
        double output[3];
        // double da = 0, db = 0, dc = 0, dd = 0, dalph = 0, dbet = 0;
        // dist_grad.execute(a, b, c, d, alph, bet, x, y, z, &da, &db, &dc, &dd, &dalph, &dbet);
        // std::cerr << "given params " << " a " << a << "b" << b << " c " << c << " d " << d << " alph " << alph << " bet " << bet << std::endl;
        double da = dist_a.execute(a, b, c, d, alph, bet, x, y, z);
        double db = dist_b.execute(a, b, c, d, alph, bet, x, y, z);
        double dc = dist_c.execute(a, b, c, d, alph, bet, x, y, z);
        double dd = dist_d.execute(a, b, c, d, alph, bet, x, y, z);
        double dalph = dist_alph.execute(a, b, c, d, alph, bet, x, y, z);
        double dbet = dist_bet.execute(a, b, c, d, alph, bet, x, y, z);
        double dist = DistanceToPoint(a, b, c, d, alph, bet, x, y, z);
        Jacobian[i * 6] = da;
        Jacobian[i * 6 + 1] = db;
        Jacobian[i * 6 + 2] = dc;
        Jacobian[i * 6 + 3] = dd;
        Jacobian[i * 6 + 4] = dalph;
        Jacobian[i * 6 + 5] = dbet;
        // std::cerr << "dist " << dist << " da " << da << " db " << db << " dc " << dc << " dd " << dd << " dalph " << dalph << " dbet " << dbet << std::endl;
    }
}

double SquareErr(double *points, int nr_of_points, double a, double b, double c, double d, double alph, double bet)
{
    // double dist[nr_of_points];
    double dist;
    double square_err = 0;
    for (int i = 0; i < nr_of_points; i++)
    {
        double x = points[i * 3];
        double y = points[i * 3 + 1];
        double z = points[i * 3 + 2];
        dist = DistanceToPoint(a, b, c, d, alph, bet, x, y, z);
        square_err += (dist * dist);
    }
    return square_err;
}

void GradientDescent(double *points, int nr_of_points)
{
    double a = 5.2122, b = 0.1, c = 0.9835, d = 1.707055, alph = -3.60384, bet = 1.13255;
    // double a = 15.7, b = 0.1, c = -34.40835, d = -1.107055, alph = -1.0384, bet = 3.3255;
    double lambda = 0.00001;
    double jacobian[nr_of_points * 6];
    double tjacobian[6 * nr_of_points];
    double dist[nr_of_points];
    double square_err = SquareErr(points, nr_of_points, a, b, c, d, alph, bet);
    double params[6] = {0};
    double prev_square_er = SquareErr(points, nr_of_points, a, b, c, d, alph, bet);
    std::cerr << square_err << std::endl;
    for (int i = 0; i < 2000; i++)
    {
        DistancesToAllPoints(points, nr_of_points, a, b, c, d, alph, bet, dist);
        Jacobian(points, nr_of_points, a, b, c, d, alph, bet, jacobian);
        Transpose(jacobian, nr_of_points, 6, tjacobian);

        double y_dist[nr_of_points];
        ScalarMultiply(dist, nr_of_points, 1, -1, y_dist);
        double h[6];
        MatrixMultiply(tjacobian, y_dist, 6, nr_of_points, 1, h);
        ScalarMultiply(h, 6, 1, lambda, h);
        // PrintMatrix("h", h, 1, 6);
        double new_square_err = SquareErr(points, nr_of_points, a + h[0], b + h[1], c + h[2], d + h[3], alph + h[4], bet + h[5]);
        // std::cerr << "lambda " << lambda << " squares distance: " << new_square_err << std::endl;

        if (new_square_err < prev_square_er)
        {
            // std::cerr << "Improvement!\n";
            lambda = lambda * 10;
        }
        else
        {
            lambda = lambda / 10;
            continue;
        }
        a += h[0];
        b += h[1];
        c += h[2];
        d += h[3];
        alph += h[4];
        bet += h[5];
        if (new_square_err < square_err)
        {
            square_err = new_square_err;
            params[0] = a;
            params[1] = b;
            params[2] = c;
            params[3] = d;
            params[4] = alph;
            params[5] = bet;
        }
        prev_square_er = new_square_err;
        std::cerr << "New params: " << a << " " << b << " " << c << " " << d << " " << alph << " " << bet << " ";
        std::cerr << "lambda: " << lambda << " squares distance: " << new_square_err << std::endl;
    }
    std::cerr << "New params: \n";
    for (int i = 0; i < 6; i++)
    {
        std::cerr << params[i] << " ";
    }
    std::cerr << square_err;
    double t = -nr_of_points / 2;
    for (int i = 0; i < 10 * nr_of_points; i++)
    {
        t += 0.1;
        double output[3];
        HelixPoint(params[0], params[1], params[2], params[3], params[4], params[5], t, output);
        double x = output[0], y = output[1], z = output[2];

        std::cout << x << " " << y << " " << z << "\n";
    }
    std::cout << "end\n";
}
//--------------------------------------------------------------------------------------------------------------------------

double Lambda(double *points, int nr_of_points, double &a, double &b, double &c, double &d, double &alph, double &bet, double lambda, double &old_square_err, double *results)
{
    // double dist[nr_of_points];
    // DistancesToAllPoints(points, nr_of_points, a + results[0], b + results[1], c + results[2], d + results[3], alph + results[4], bet + results[5], dist);
    // double square_err = 0;
    double new_lambda;
    // for (int i = 0; i < nr_of_points; i++)
    // {
    //     square_err += (dist[i] * dist[i]);
    // }
    double square_err = SquareErr(points, nr_of_points, a + results[0], b + results[1], c + results[2], d + results[3], alph + results[4], bet + results[5]);
    std::cerr << "SQUARE ERR " << square_err << std::endl;
    if ((square_err >= old_square_err) && (lambda < 1000))
        new_lambda = lambda * 10;
    else
    {
        std::cerr << "IMPROVEMENTS!";
        a += results[0];
        b += results[1];
        c += results[2];
        d += results[3];
        alph += results[4];
        bet += results[5];
        new_lambda = lambda / 10;
        old_square_err = square_err;
    }
    return new_lambda;
}
void LevenbergMarquardt(double *points, int nr_of_points)
{

    auto seed = time(nullptr);
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> uniform(-2 * std::numbers::pi_v<double>, 2 * std::numbers::pi_v<double>);
    // double a = uniform(rng), b = uniform(rng), c = uniform(rng), d = uniform(rng), alph = uniform(rng), bet = uniform(rng);
    // double a = 5.122, b = -4.79395, c = -26.0835, d = -4.07055, alph = -3.60384, bet = 1.13255;
    // double a = 5.2122, b = 2, c = 0.40835, d = 0.707055, alph = -3.60384, bet = 1.13255;
    // double a = 20.5, b = 0.001, c = 40.0, d = 15.207055, alph = -3.60384, bet = 1.13255;
    // double a = 6.2122, b = 0.1, c = 15.9835, d = 11.707055, alph = -3.60384, bet = 1.13255;
    double a = 6.2122, b = 0.1, c = 1.9835, d = 1.707055, alph = -3.60384, bet = 1.13255;

    int diff_params = 6;
    double lambda = 1;
    double lambda_change = 1;
    double old_square_err;
    double jacobian[nr_of_points * diff_params];
    double tjacobian[diff_params * nr_of_points];
    double tjj[diff_params * diff_params];
    double results[diff_params];
    {
        double dist[nr_of_points];
        DistancesToAllPoints(points, nr_of_points, a, b, c, d, alph, bet, dist);
        old_square_err = 0;
        for (int i = 0; i < nr_of_points; i++)
        {
            old_square_err += (dist[i] * dist[i]);
        }
    }

    //  --------------------------------------------------------------------------------------------------------------------
    for (int i = 0; i < 1000; i++)
    {
        // if (old_square_err < 0.1)
        //     break;

        Jacobian(points, nr_of_points, a, b, c, d, alph, bet, jacobian);
        // PrintMatrix("jacobian", jacobian, nr_of_points, 5);

        Transpose(jacobian, nr_of_points, diff_params, tjacobian);

        MatrixMultiply(tjacobian, jacobian, diff_params, nr_of_points, diff_params, tjj);

        double diag[diff_params * diff_params];
        DiagOfSquareM(tjj, diff_params, diag);

        double identity[diff_params * diff_params];
        ScalarMultiply(diag, diff_params, diff_params, lambda, identity);
        double left_side[diff_params * diff_params];
        AddMatrices(tjj, identity, diff_params, diff_params, left_side);
        double dist[nr_of_points];
        DistancesToAllPoints(points, nr_of_points, a, b, c, d, alph, bet, dist);
        double right_side[diff_params * 1];
        MatrixMultiply(tjacobian, dist, diff_params, nr_of_points, 1, right_side);
        ScalarMultiply(right_side, 1, diff_params, -1, right_side);

        // left side is 6x6, right side is 6x1, so h is 6x1.
        double forward_elim[diff_params * diff_params];
        double unchanged_rs[diff_params];
        CopyMatrix(right_side, diff_params, unchanged_rs);
        PrintMatrix("left_side", left_side, diff_params, diff_params);
        PrintMatrix("right_side", unchanged_rs, 1, diff_params);
        ForwardElim(left_side, diff_params, right_side, forward_elim);
        PrintMatrix("right_side", right_side, 1, diff_params);
        BackSub(forward_elim, diff_params, right_side, results);
        CheckSolution(left_side, diff_params, unchanged_rs, results);
        lambda = Lambda(points, nr_of_points, a, b, c, d, alph, bet, lambda, old_square_err, results);
        std::cerr << "New params: " << a << " " << b << " " << c << " " << d << " " << alph << " " << bet << " ";
        std::cerr << "lambda: " << lambda << " squares distance: " << old_square_err << std::endl;
    }
    double t = -nr_of_points / 2;
    for (int i = 0; i < 10 * nr_of_points; i++)
    {
        t += 0.1;
        double output[3];
        HelixPoint(a, b, c, d, alph, bet, t, output);
        double x = output[0], y = output[1], z = output[2];

        std::cout << x << " " << y << " " << z << "\n";
    }
    std::cout << "end\n";
}
