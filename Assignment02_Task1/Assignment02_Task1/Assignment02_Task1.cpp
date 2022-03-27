#include <iostream>
#include <omp.h>
#include <cstdio>
#include <stdio.h>
#include <string.h>
#include <cstdlib>
#include <cstdio>
#define _USE_MATH_DEFINES
#include <math.h>

void AllocateVector(double** Vector, int size)
{
    (*Vector) = new double[size];
}

void FreeVector(double** Vector)
{
    delete[](*Vector);
}

void PrintVector(double* Vector, int size)
{
    for (int i = 0; i < size; i++) {
        printf("%f", Vector[i]);
        printf("\n");
    }
}

double mu(double t) { // at x = 0 or x = N
    return 0;
}

double u0(double x) {  // at t = 0 
    return 0.1 * sin(M_PI * x);
}

double alpha(double x) {
    return exp(cos(x))/ 10;
}

void explicit_solu(double* solution, double tau, double h, int L, int N, double T) {
    double x;
    double t;
    for (int j = 0; j < L; j++) {
        t = j * tau;
        solution[j * N] = mu(t);
        solution[((j + 1) * N) - 1] = mu(t);
    }
    for (int i = 1; i < N - 1; i++) {
        x = i * h;
        solution[i] = u0(x);
        solution[N + i] = u0(x);
    }
    for (int j = 1; j < L - 1; j++) {

        for (int i = 1; i < N - 1; i++) {

            solution[(j + 1) * N + i] = (2 * solution[N * j + i]) - solution[N * (j - 1) + i] + (pow(tau, 2) * alpha(x) * (1 / pow(h, 2)) * (solution[(N * j) + i + 1] - 2 * solution[(N * j) + i] + solution[(N * j )+ i - 1]));
        }
    }
    PrintVector(solution, N * L);
}

int main()
{
    int L = 20;
    int N = 100;
    double T = 2;
    double tau = (double)T/L;
    double h = (double)1/N;
    double *solu;
    AllocateVector(&solu, N * L);
    explicit_solu(solu, tau, h, L, N, T);
    return 0;
}