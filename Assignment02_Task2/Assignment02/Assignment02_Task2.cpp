// Task 2 
#include "math.h"
#include <omp.h>
#include <cstdio>
#include <stdio.h>
#include <string.h>
#include <cstdlib>
#include <cstdio>
using namespace std;
#define _USE_MATH_DEFINES
#define EPSILON 0.001
#define LEFT_BOUND 1
#define RIGHT_BOUND 1
#define N_MAX 100000

double f(double x, double y) {
	return -4;
}
double mu1(double y) {
	return (1 + pow(y+2, 2));
}
double mu2(double y) {
	return (4 + pow(y+2, 2));
}
double mu3(double x) {
	return (pow(x+1, 2) + 4);
}
double mu4(double x) {
	return (pow(x+1,2) + 9);
}
double FuncU(double x, double y) {
	return (pow(x+1,2) + pow(y+2,2));
}

// Vector allocation
void AllocateVector(double** Vector, int size)
{
	(*Vector) = new double[size];
}
// Vector release
void FreeVector(double** Vector)
{
	delete[](*Vector);
}
// Vector printing
void PrintVector(char* VectorName, double* Vector, int size)
{
	printf("%s\n", VectorName);
	for (int i = 0; i < size; i++)
		printf("%lf ", Vector[i]);
	printf("\n");
}
//--------------------------------------------------------------------------
// Initialization of the data by the partial differential equation
//--------------------------------------------------------------------------
// matrix allocation (n, m)
// Diagonals of the matrix are stored by rows
// the shift if the diagonal elements is stored in the index array
void CreateDUMatrix(int n, int m, double** Matrix, double** Index)
{
	double hsqr = (double)n * n / LEFT_BOUND / LEFT_BOUND; // 1/h
	double ksqr = (double)m * m / RIGHT_BOUND / RIGHT_BOUND; // 1/k
	double A = 2 * (hsqr + ksqr);
	int size = (n - 1) * (m - 1), bandWidth = 5;
	AllocateVector(Matrix, size * bandWidth);
	AllocateVector(Index, bandWidth);
	(*Index)[0] = -n + 1; (*Index)[1] = -1; (*Index)[2] = 0; (*Index)[3] = 1; (*Index)[4] = n - 1;
	for (int i = 0; i < size; i++)
	{
		if (i >= n - 1) (*Matrix)[i * bandWidth] = -ksqr;
		else (*Matrix)[i * bandWidth] = 0.0;
		if (i % (n - 1) != 0) (*Matrix)[i * bandWidth + 1] = -hsqr;	
		else (*Matrix)[i * bandWidth + 1] = 0.0;
		(*Matrix)[i * bandWidth + 2] = A;
		if ((i + 1) % (n - 1) != 0) (*Matrix)[i * bandWidth + 3] = -hsqr;
		else (*Matrix)[i * bandWidth + 3] = 0.0;
		if (i < (n - 1) * (m - 2)) (*Matrix)[i * bandWidth + 4] = -ksqr;
		else (*Matrix)[i * bandWidth + 4] = 0.0;
	}
}
// The right-hand side of the equation
void CreateDUVector(int n, int m, double** Vector)
{
	double h = LEFT_BOUND / (double)n;
	double k = RIGHT_BOUND / (double)m;
	double hsqr = (double)n * n / LEFT_BOUND / LEFT_BOUND;
	double ksqr = (double)m * m / RIGHT_BOUND / RIGHT_BOUND;
	AllocateVector(Vector, (n - 1) * (m - 1));
	for (int j = 0; j < m - 1; j++)
	{
		for (int i = 0; i < n - 1; i++)
			(*Vector)[j * (n - 1) + i] = f((double)(i) * h, (double)(j) * k);

		(*Vector)[j * (n - 1)] += hsqr * mu1((double)(j) * k);
		(*Vector)[j * (n - 1) + n - 2] += hsqr * mu2((double)(j) * k);
	}
	for (int i = 0; i < n - 1; i++)
	{
		(*Vector)[i] += ksqr * mu3((double)(i) * h);
		(*Vector)[(m - 2) * (n - 1) + i] += ksqr * mu4((double)(i) * h);
	}
}
void GetFirstApproximation(double* Result, int size)
{
	for (int i = 0; i < size; i++)
		Result[i] = 0.0;
}

double BandOverRelaxation(double* Matrix, double* Vector, double* Result, double* Index, int
	size, int bandWidth,
	double WParam, double Accuracy, int& StepCount)
{
	double CurrError; 
	double sum, TempError;
	int k;
	int ii, index = Index[bandWidth - 1], bandHalf = (bandWidth - 1) / 2;
	StepCount = 0;
	do
	{
		CurrError = -1.0;
		for (int i = index; i < size + index; i++)
		{
			ii = i - index;
			TempError = Result[i];
			sum = 0.0;   
			for (int j = 0; j < bandWidth; j++) {
				k = i + Index[j];
				sum += Matrix[(ii * bandWidth) + j] * Result[k];
			}
			Result[i] = (Vector[ii] - sum) * WParam / Matrix[ii * bandWidth + bandHalf] + Result[i];
			TempError = fabs(Result[i] - TempError);
			if (TempError > CurrError) CurrError = TempError;
		}
		StepCount++;
	} while ((CurrError > Accuracy) && (StepCount < N_MAX));
	return CurrError;
}
double SolvePoisson(int n, int m, double* Solution, double Accuracy, double& ORAccuracy,
	int& StepCount)
{
	double* Matrix, * Vector, * Result;
	double* Index;
	int size = (n - 1) * (m - 1), ResSize = size + 2 * (n - 1), bandWidth = 5;
	double start, finish; double time;
	double WParam, step = (n / LEFT_BOUND > m / RIGHT_BOUND) ?
		(double)LEFT_BOUND / n : (double)RIGHT_BOUND / m;
	CreateDUMatrix(n, m, &Matrix, &Index);
	CreateDUVector(n, m, &Vector);
	AllocateVector(&Result, ResSize);
	GetFirstApproximation(Result, ResSize);
	WParam = 1.5;
	start = omp_get_wtime();
	ORAccuracy = BandOverRelaxation(Matrix, Vector, Result, Index, size, bandWidth,
		WParam, Accuracy, StepCount);
	finish = omp_get_wtime();
	time = (finish - start);
	memcpy(Solution, Result + n - 1, sizeof(double) * size);
	FreeVector(&Matrix);
	FreeVector(&Index);
	FreeVector(&Vector);
	FreeVector(&Result);
	return time;
}
double SolutionCheck(double* solution, int n, int m)
{
	double h = LEFT_BOUND / (double)n, k = RIGHT_BOUND / (double)m;
	double err = 0, temp;
	for (int j = 0; j < m - 1; j++)
		for (int i = 0; i < n - 1; i++)
		{
			temp = fabs(solution[j * (n - 1) + i] - FuncU((double)(i) * h, (double)(j) * k));
			if (temp > err)
				err = temp;
		}
	return err;
}
void plotSolution(int n) {
	double h = LEFT_BOUND / (double)n;
	double k = RIGHT_BOUND / (double)n;
	FILE* gnuplot = _popen("C:\Program Files\gnuplot\bin\gnuplot -persist" ,"w");
	fprintf(gnuplot, "set view map\n");
	fprintf(gnuplot, "set dgrid3d\n");
	fprintf(gnuplot, "set pm3d interpolate 0,0\n");
	fprintf(gnuplot, "set title 'Heatmap'\n");
	fprintf(gnuplot, "set xrange [-0.001:1.001]\n");
	fprintf(gnuplot, "splot '-' using 1:2:3 with pm3d\n");
	for (int i = 0; i <= n; ++i) {
		double x = h * i;
		for (int j = 0; j <= n; ++j) {
			double y = k * j;
			fprintf(gnuplot, "%f %f %f\n", x, y + 2, FuncU(x, y));
		}
	}
	fprintf(gnuplot, "e\n");
	_pclose(gnuplot);
}
int main(int argc, char* argv[])
{
	int n[5] = { 10, 50, 100, 500, 1000 };
	int m[5] = { 10, 50, 100, 500, 1000 };
	double Epsilon[3] = { 0.001, 0.0001, 0.00001 };
	printf("Nabila Adawy, n.roshdy@innopolis.university, PA-2 (problem 2) \n");
	for (int i = 0; i < 5; i++) {
		printf("\n\n|  Grid Size  |  Estimated Accuracy  |  Iter Number  |   Achieved Accuracy   |\n");
	
		for (int j = 0; j < 3; j++) {
			int StepCount;
			int size;
			double time;
			double Accuracy = Epsilon[j];
			double AcAccuracy;
			double Correctness;
			double* Solution;
			size = (n[i] - 1) * (m[i] - 1);
			AllocateVector(&Solution, size);
			time = SolvePoisson(n[i], m[i], Solution, Accuracy, AcAccuracy, StepCount);
			Correctness = SolutionCheck(Solution, n[i], m[i]);
			printf("|    %5i    |       %5f       |   %7i     |  %15f      |\n", n[i], Epsilon[j], StepCount, AcAccuracy);
			FreeVector(&Solution);
		}
	}
	plotSolution(n[2]);
	return 0;
}