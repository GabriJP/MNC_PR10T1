#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <mkl.h>
#include "MDF.h"

void getparams(const int N, const double* x, double* params)
{
	auto A = params;
	auto B = &(params[N]);
	auto C = &(params[2 * N]);
	auto D = &(params[3 * N]);
	for (auto i = 0; i < N; ++i)
	{
		A[i] = x[i] * x[i];
		B[i] = x[i];
		C[i] = -1.0;
		D[i] = 3.0 * x[i] * x[i];
	}
}

double exactSolution(double x)
{
	return x + 1.0 / x + x * x;
}

int main(int arc, char* argv[])
{
	auto mdf = MDF();
	auto a = 0.5;
	auto b = 2.0;
	auto N = 150;
	auto h = (b - a) / static_cast<double>(N - 1);
	auto x = static_cast<double*>(mkl_malloc(N * sizeof(double), 64));
	auto y = static_cast<double*>(mkl_malloc(N * sizeof(double), 64));
	auto params = static_cast<double*>(mkl_malloc(4 * N * sizeof(double), 64));
	for (auto i = 0; i < N; ++i)
	{
		x[i] = a + static_cast<double>(i * h);
	}
	getparams(N, x, params);
	mdf.solve(N, h, y, params, CONTOUR_TYPE1, exactSolution(a), exactSolution(b));
	auto yexact = static_cast<double*>(mkl_malloc(N * sizeof(double), 64));
	for (auto i = 0; i < N; ++i)
	{
		yexact[i] = exactSolution(x[i]);
	}
	auto error = 0.0;
	for (auto i = 0; i < N; ++i)
	{
		error += fabs((yexact[i] - y[i]) / yexact[i]);
	}
	error /= N;
	printf("\nError relativo promedio(%%): %g\n", 100.0 * error);
	mkl_free(x);
	mkl_free(y);
	mkl_free(yexact);
	getchar();
	return 0;
}