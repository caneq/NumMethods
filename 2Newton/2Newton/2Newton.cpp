#include "pch.h"
#include <iostream>
#include <iomanip>
using namespace std;


int gaus(double A[][2], double* B, int n, double* res);
double f1(double x1, double x2) {
	return 2 * x1*x1*x1 - x2 * x2 - 1;
}

double f2(double x1, double x2) {
	return x1 * x2*x2*x2 - x2 - 4;
}

double df1x1(double x1, double x2) {
	return 6 * x1*x1;
}
double df1x2(double x1, double x2) {
	return -2 * x2;
}

double df2x1(double x1, double x2) {
	return x2 * x2*x2;
}
double df2x2(double x1, double x2) {
	return 3 * x1 * x2*x2 - 1;
}


void calcJAn(double J[][2], const double &x1, const double &x2) {
	J[0][0] = df1x1(x1, x2);
	J[0][1] = df1x2(x1, x2);
	J[1][0] = df2x1(x1, x2);
	J[1][1] = df2x2(x1, x2);
}

void calcJInc(double J[][2], const double &x1, const double &x2) {
	const double M = 0.01;
	const double x1M = x1 * M;
	const double x2M = x2 * M;
	J[0][0] = (f1(x1 + x1M, x2) - f1(x1, x2)) / x1M;
	J[0][1] = (f1(x1, x2 + x2M) - f1(x1, x2)) / x2M;
	J[1][0] = (f2(x1 + x1M, x2) - f2(x1, x2)) / x1M;
	J[1][1] = (f2(x1, x2 + x2M) - f2(x1, x2)) / x2M;
}

int newton(double& x1, double& x2, double e1, double e2, int NIT) {
	int k = 1;
	cout << left << setw(4) << "k" << setw(20) << "d1" << setw(20) << "d2" << setw(20) << "x1" << setw(20) << "x2" << endl;
	double F[2], J[2][2];
	double dx[2] = { 0 };
	double x1k, x2k;
	double d1, d2;
	double tmp;
	do {
		F[0] = -f1(x1, x2); F[1] = -f2(x1, x2);
		calcJInc(J, x1, x2);

		gaus(J, F, 2, dx);


		x1k = x1 + dx[0];
		x2k = x2 + dx[1];

		d1 = fabs(f1(x1, x2));
		tmp = fabs(f2(x1, x2));
		if (tmp > d1) d1 = tmp;

		d2 = fabs(x1k - x1) / (x1k >= 1 ? x1k : 1);
		tmp = fabs(x2k - x2) / (x2k >= 1 ? x2k : 1);
		if (tmp > d2) d2 = tmp;

		x1 = x1k;
		x2 = x2k;

		cout << left << setw(4) << k << setw(20) << d1 << setw(20) << d2 << setw(20) << x1 << setw(20) << x2 << endl;
		if (k >= NIT) {
			cout << "IER = 2\n";
			return 2;
		}
		k++;
	} while (d1 > e1 && d2 > e2);

	return 0;
}

int main() {
	double x1 = 1;
	double x2 = 1;
	newton(x1, x2, 1e-9, 1e-9, 1000);
	cout << "\nx1 = " << setprecision(15) << setw(26) << x1 << "x2 = " << setw(26) << x2 << endl;
}

int gaus(double A[][2], double* B, int n, double* res) {

	double* B1 = new double[n];
	for (int i = 0; i < n; i++) {
		B1[i] = B[i];
	}
	double** A1 = new double*[n];
	for (int i = 0; i < n; i++) {
		A1[i] = new double[n];
	}
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			A1[i][j] = A[i][j];
		}
	}

	//Прямой ход:
	for (int k = 0; k < n; k++) {
		int max = k;
		for (int i = k + 1; i < n; i++) {
			if (fabs(A1[i][k]) > fabs(A1[max][k])) max = i;
		}
		if (k != max) {
			swap(A1[k], A1[max]);
			swap(B1[k], B1[max]);
			//print(A1, B1, n);
		}
		double AMAIN = A1[k][k];
		if (fabs(AMAIN) <= DBL_EPSILON) {
			for (int i = 0; i < n; i++) delete[] A1[i];
			delete[] A1;
			delete[] B1;
			return 1; //корней нет или их бесконечно много
		}

		for (int j = k; j < n; j++) {
			A1[k][j] /= AMAIN;
		}
		B1[k] /= AMAIN;
		//print(A1, B1, n);
		for (int i = k + 1; i < n; i++) {
			double kef = A1[i][k];
			for (int j = k; j < n; j++) {
				A1[i][j] -= A1[k][j] * kef;
			}
			B1[i] -= B1[k] * kef;
		}
		//print(A1, B1, n);
	}

	//Обратный ход:
	for (int i = n - 1; i >= 0; i--) {
		res[i] = B1[i];
		for (int j = i + 1; j < n; j++) {
			res[i] -= A1[i][j] * res[j];
		}
	}
	for (int i = 0; i < n; i++) delete[] A1[i];
	delete[] A1;
	delete[] B1;
	return 0;
}