#pragma once
#include <cassert>
#include <iostream>

// Matrices are written as 1D arrays!
inline void MatrixMultiply(double *a, double *b, int arows, int acols, int bcols, double *output)
{
    for (int i = 0; i < arows; i++)
    {
        for (int j = 0; j < bcols; j++)
        {
            double sum = 0;
            for (int k = 0; k < acols; k++)
                sum = sum + a[i * acols + k] * b[k * bcols + j];
            output[i * bcols + j] = sum;
        }
    }
}

inline void Transpose(double *input, int rows, int cols, double *output)
{
    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            int i_input = i * cols + j;

            int i_output = j * rows + i;

            output[i_output] = input[i_input];
        }
    }
}

inline void DiagOfSquareM(double *input, int height, double *diag)
{
    // Works for square matrices only
    for (int i = 0; i < height * height; i++)
    {
        diag[i] = 0;
    }
    for (int i = 0; i < height; i++)
    {
        diag[i * height + i] = input[i * height + i];
    }
}

inline void ScalarMultiply(double *matrix, int rows, int cols, double number, double *output)
{
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            output[i * cols + j] = number * matrix[i * cols + j];
        }
    }
}
inline void CopyMatrix(double *matrix, int size, double *output)
{
    for (int i = 0; i < size; i++)
    {
        output[i] = matrix[i];
    }
}
inline void AddMatrices(double *a, double *b, int rows, int cols, double *output)
{
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            output[i * cols + j] = a[i * cols + j] + b[i * cols + j];
        }
    }
}

inline void swap_row(double matrix[6][6], int i, int j)
{

    for (int k = 0; k < 6; k++)
    {
        double temp = matrix[i][k];
        matrix[i][k] = matrix[j][k];
        matrix[j][k] = temp;
    }
}
inline void ForwardElim(double input[6][6], double *res, double output[6][6])
{
    for (int i = 0; i < 6; ++i)
    {
        for (int j = 0; j < 6; ++j)
        {
            output[i][j] = input[i][j];
        }
    }
    for (int i = 0; i < 6; i++)
    {
        int i_max = i;
        double v_max = output[i_max][i];

        for (int j = i + 1; j < 6; j++)
            if (std::abs(output[j][i]) > std::abs(v_max) && output[j][i] != 0)
                v_max = output[j][i], i_max = j;
        if (i_max != i)
        {
            swap_row(output, i, i_max);
            double temp = res[i];
            res[i] = res[i_max];
            res[i_max] = temp;
        }

        if (output[i][i] == 0.0)
        {
            std::cerr << "Mathematical Error!";
            std::cerr << "Input that caused the error is:";
            for (int i = 0; i < 6; i++)
            {
                for (int j = 0; j < 6; j++)
                {
                    std::cerr << output[i][j] << " ";
                }
                std::cerr << std::endl;
            }
        }
        for (int j = i + 1; j < 6; j++)
        {
            double ratio = output[j][i] / output[i][i];

            for (int k = 0; k < 6; k++)
            {
                output[j][k] = output[j][k] - ratio * output[i][k];
                if (std::abs(output[j][k]) <= 1e-15)
                {
                    output[j][k] = 0;
                }
            }
            res[j] = res[j] - ratio * res[i];
            if (std::abs(res[j]) <= 1e-15)
            {
                res[j] = 0;
            }
        }
    }

    // TO DO. Take this out
    std::cerr << "Forward elimination result:" << std::endl;
    for (int j = 0; j < 6; j++)
    {
        for (int k = 0; k < 6; k++)
        {
            std::cerr << output[j][k] << " ";
        }
        std::cerr << std::endl;
    }
}

void BackSub(double input[6][6], double *right_side, double *results)
{
    /*Back substitution and the result of Gaussian elimination*/
    for (int i = 5; i > -1; i--)
    {
        results[i] = right_side[i];
        for (int j = 5; j > i; j--)
        {
            results[i] -= input[i][j] * results[j];
        }
        results[i] /= input[i][i];
        if (std::abs(results[i]) <= 1e-15)
        {
            results[i] = 0;
        }
    }
}
void CheckSolution(double input[6][6], double *right_side, double *results)
{
    for (int i = 0; i < 6; i++)
    {
        double sum = 0;
        for (int j = 0; j < 6; j++)
        {
            sum += input[i][j] * results[j];
        }
        if (std::abs(sum - right_side[i]) >= 1e-5)
        {
            std::cerr << "Wrong solution " << sum << " " << right_side[i] << std::endl;
            std::abort();
        }
    }
}