/*
Método de las Diferencias Finitas
ULPGC, EII, MNC
*/

#include <cstdio>
#include <cstdlib>

#include <mkl.h>

#include "MDF.h"

MDF::MDF(){}

MDF::~MDF(){}

int MDF::solve(const int N, const double h,
	double *y, const double *params, const int type, 
	const double ca, const double cb){

	double *RHS = (double *)mkl_malloc(N*sizeof(double), 64);
	double *LHS = (double *)mkl_malloc(3 * N*sizeof(double), 64);
	int *ipiv = (int *)mkl_malloc(N*sizeof(int),32);

	ensamble_kernel(N, h, params, LHS, RHS);
	switch (type){
	case CONTOUR_TYPE1:
		ensamble_contour1(N, h, LHS, RHS, ca, cb);
		break;
	case CONTOUR_TYPE2:
		ensamble_contour2(N, h, LHS, RHS, ca, cb);
		break;
	case CONTOUR_TYPE3:
		ensamble_contour3(N, h, LHS, RHS, ca, cb);
		break;
	case CONTOUR_TYPE4:
		ensamble_contour4(N, h, LHS, RHS, ca, cb);
		break;
	default:
		break;
	}

	int info = LAPACKE_dgtsv(LAPACK_ROW_MAJOR, N, 1, &(LHS[2*N]), &(LHS[N]), &(LHS[1]), RHS, 1);
	
	for (int i = 0; i < N; i++) y[i] = RHS[i];

	mkl_free(LHS);
	mkl_free(RHS);
	mkl_free(ipiv);

	return info;
}

// Ensambla todo menos las condiciones de contorno
void MDF::ensamble_kernel(const int N, const double h, const double *params,
	double *LHS, double *RHS){

	const double *A = params;
	const double *B = &(params[N]);
	const double *C = &(params[2*N]);
	const double *D = &(params[3*N]);

	double *E, *F, *G, *H;
	E = (double *)mkl_malloc(N*sizeof(double), 64);
	F = (double *)mkl_malloc(N*sizeof(double), 64);
	G = (double *)mkl_malloc(N*sizeof(double), 64);
	H = (double *)mkl_malloc(N*sizeof(double), 64);

	for (int i = 1; i < N-1; i++){
		E[i] = A[i];
		F[i] = h*h*C[i] - h*B[i] - 2*A[i];
		G[i] = A[i] + h*B[i];
		H[i] = h*h*D[i];
	}

	LHS[0] = LHS[1] = LHS[N] = 0.0;
	for (int i = 2; i < N; i++)     LHS[i] = G[i - 1];
	for (int i = 1; i < N - 1; i++) LHS[N+i] = F[i];
	for (int i = 0; i < N - 2; i++) LHS[2*N+i] = E[i + 1];
	LHS[2*N-1] = LHS[3*N-1] = LHS[3*N-2] = 0.0;

	for (int i = 1; i< N - 1; i++) RHS[i] = H[i];
	RHS[0] = RHS[N - 1] = 0.0;

	mkl_free(E);
	mkl_free(F);
	mkl_free(G);
	mkl_free(H);
}

void MDF::ensamble_contour1(const int N, const double h, double *LHS, double *RHS, double ya, double yb){

	// condiciones en a
	LHS[N] =  1.0;
	RHS[0] = ya;

	// condiciones en b
	LHS[2 * N - 1] = 1.0;
	RHS[N - 1] = yb;
}

void MDF::ensamble_contour2(const int N, const double h, double *LHS, double *RHS, double ya, double y1b){

	//condiciones en a
	LHS[N] = 1.0;
	RHS[0] = ya;

	// condiciones en b
	LHS[2 * N - 1] = 1.0;
	LHS[3 * N - 2] = -1.0;
	RHS[N - 1] = y1b*h;
}

void MDF::ensamble_contour3(const int N, const double h, double *LHS, double *RHS, double y1a, double yb){

	//condiciones en a
	LHS[N] = 1.0;
	LHS[2] = -1.0;
	RHS[0] = y1a*h;

	// condiciones en b
	LHS[2 * N - 1] = 1.0;
	RHS[N - 1] = yb;
}

void MDF::ensamble_contour4(const int N, const double h, double *LHS, double *RHS, double y1a, double y1b){

	//condiciones en a
	LHS[N] = 1.0;
	LHS[2] = -1.0;
	RHS[0] = y1a*h;

	// condiciones en b
	LHS[2 * N - 1] = 1.0;
	LHS[3 * N - 2] = -1.0;
	RHS[N - 1] = y1b*h;
}
