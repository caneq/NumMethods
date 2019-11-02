#include "pch.h"
#include <iostream>
#include <vector>
#include <iomanip>
#include <string>
using namespace std;

typedef double(*func)(const vector<double>&);

double f1(const vector<double>& x) {
	return 2 * x[0] * x[0] * x[0] - x[1] * x[1] - 1;
}

double f2(const vector<double>& x) {
	return x[0] * x[1] * x[1] * x[1] - x[1] - 4;
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

int gaus(vector<vector<double>> A, vector<double> B, vector<double>& res) {
	int n = B.size();
	if (res.size() != B.size()) {
		res = B;
	}

	for (int k = 0; k < n; k++) {
		int max = k;
		for (int i = k + 1; i < n; i++) {
			if (fabs(A[i][k]) > fabs(A[max][k])) max = i;
		}
		if (k != max) {
			swap(A[k], A[max]);
			swap(B[k], B[max]);
		}
		double AMAIN = A[k][k];
		if (fabs(AMAIN) <= DBL_EPSILON) {
			return 1;
		}

		for (int j = k; j < n; j++) {
			A[k][j] /= AMAIN;
		}
		B[k] /= AMAIN;
		for (int i = k + 1; i < n; i++) {
			double kef = A[i][k];
			for (int j = k; j < n; j++) {
				A[i][j] -= A[k][j] * kef;
			}
			B[i] -= B[k] * kef;
		}
	}
	for (int i = n - 1; i >= 0; i--) {
		res[i] = B[i];
		for (int j = i + 1; j < n; j++) {
			res[i] -= A[i][j] * res[j];
		}
	}
	return 0;
}

void calcJAn(double J[][2], const double &x1, const double &x2) {
	J[0][0] = df1x1(x1, x2);
	J[0][1] = df1x2(x1, x2);
	J[1][0] = df2x1(x1, x2);
	J[1][1] = df2x2(x1, x2);
}

void calcJInc(vector<vector<double>>& J, const vector<func> &f, vector<double> x) {
	J.clear();
	J.reserve(f.size());
	const double M = 0.01;
	vector<double> tmp;
	for (int i = 0; i < f.size(); i++) J.push_back(tmp);
	for (int i = 0; i < f.size(); i++) {
		double xM = x[i] * M;
		for (int j = 0; j < f.size(); j++) {
			double fx = f[j](x);
			x[i] += xM;
			J[j].push_back((f[j](x) - fx) / xM);
			x[i] -= xM;
		}
	}
}

int newton(const vector<func> &f, vector<double>& x, double e1, double e2, int NIT) {
	cout << left << setw(4) << "k" << setw(20) << "d1" << setw(20) << "d2";
	string xString("x");
	for (int i = 0; i < f.size(); i++) {
		cout << setw(20) << (xString + to_string(i));
	}
	cout << endl;
	int k = 1;
	vector<double> F = x;
	if (x.size() != f.size()) {
		x.clear();
		x.reserve(f.size());
		for (int i = 0; i < f.size(); i++) {
			x.push_back(0.0);
		}
	}
	vector<vector<double>> J;
	vector<double> dx;
	vector<double> resk = x;
	double d1, d2;
	do {
		for (int i = 0; i < f.size(); i++) {
			F[i] = -f[i](x);
		}
		calcJInc(J, f, x);

		gaus(J, F, dx);

		for (int i = 0; i < resk.size(); i++) {
			resk[i] = x[i] + dx[i];
		}

		d1 = DBL_MIN;
		for (int i = 0; i < x.size(); i++) {
			double tmp = fabs(f[i](x));
			if (tmp > d1) d1 = tmp;
		}

		d2 = DBL_MIN;
		for (int i = 0; i < resk.size(); i++) {
			double tmp = fabs(resk[i] - x[i]);
			if (resk[i] >= 1) {
				tmp /= resk[i];
			}

			if (tmp > d2) d2 = tmp;
		}

		x = resk;

		cout << left << setw(4) << k << setw(20) << d1 << setw(20) << d2;
		for (int i = 0; i < f.size(); i++) {
			cout << setw(20) << x[i];
		}
		cout << endl;
		if (k >= NIT) {
			cout << "IER = 2\n";
			return 2;
		}
		k++;
	} while (d1 > e1 && d2 > e2);

	return 0;
}

int main() {
	vector<func> f; f.push_back(f1); f.push_back(f2);
	vector<double> x; x.push_back(1); x.push_back(1);
	newton(f, x, 1e-10, 1e-10, 1000);
	cout << "\nx1 = " << setprecision(15) << setw(26) << x[0] << "x2 = " << setw(26) << x[1] << endl;
	return 0;
}