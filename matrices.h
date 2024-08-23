#pragma once
// To do. Remove some constants and use variables instead.

void MatrixMultiply(double a[6][200], double b[200][6], double **output)
{
    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < 6; j++)
        {
            output[i][j] = 0.0;
            for (int k = 0; k < 200; k++)
            {
                output[i][j] += a[i][k] * b[k][j];
            }
        }
    }
}

void Transpose(double input[200][6], double output[6][200])
{
    for (int i = 0; i < 200; i++)
    {
        for (int j = 0; j < 6; j++)
        {
            output[j][i] = input[i][j];
        }
    }
}

void Diag(double input[6][6], double diag[6][6])
{
    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < 6; j++)
        {
            if (i == j)
                diag[i][j] = input[i][i];
            else
                diag[i][j] = 0;
        }
    }
}

void ScalarMultiply(double matrix[6][6], double nr, double output[6][6])
{
    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < 6; j++)
        {
            output[i][j] = matrix[i][j] * nr;
        }
    }
}