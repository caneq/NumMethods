#pragma once
#include <iostream>
#include <vector>
#include <iomanip>
using namespace std;

typedef double(*func)(const vector<double> &, double);

void explicitEulerMethod(const vector<func> &f, double T, double epsI,
	double tauMax, const vector<double> &u, vector<vector<double>> &points) {
	double tk = 0;
	vector<double> yk = u;
	int N = yk.size();
	for (int i = 0; i < N; i++) {
		cout << setw(13) << "y" << i + 1;
	}
	cout << setw(13) << "tk\n";
	while (tk < T) {
		double tauK = DBL_MAX;
		for (int i = 0; i < N; i++) {
			double tmp = epsI / (fabs(f[i](yk, tk)) + epsI / tauMax);
			if (tmp < tauK) tauK = tmp;
		}

		for (int i = 0; i < N; i++) yk[i] += tauK * f[i](yk, tk);
		tk += tauK;

		vector<double> point;
		point.push_back(tk);
		for (int i = 0; i < N; i++) {
			point.push_back(yk[i]);
			//cout << setw(14) << yk[i];
		}
		points.push_back(point);
		//cout << setw(14) << tk << endl;
	}
}

int newton(const vector<func> &f, const vector<double>& yk, vector<double> &ykPlus1,
	double tauK, double tKplus1, double e1, double e2, int NIT);

double calcTauKplus1Opt(const double& epsI, vector<double> epsK, const double& tauK) {
	double res = DBL_MAX;
	for (int i = 0; i < epsK.size(); i++) {
		double tmp = sqrt(epsI / fabs(epsK[i])) * tauK;
		if (tmp < res) res = tmp;
	}
	return res;
}

double calcTauKplus1Three(const double& epsI, vector<double> epsK, const double& tauK) {
	double res = DBL_MAX;
	const double epsIDiv4 = epsI / 4;
	for (int i = 0; i < epsK.size(); i++) {
		double tmp;
		double absEpsKi = fabs(epsK[i]);
		if (absEpsKi > epsI) {
			tmp = tauK / 2;
		}
		else if (epsIDiv4 < absEpsKi && absEpsKi <= epsI) {
			tmp = tauK;
		}
		else {
			tmp = 2 * tauK;
		}
		if (tmp < res) res = tmp;
	}

	return res;
}

void implicitEulerMethod(const vector<func> &f, double T, double epsI, double tauMin,
	double tauMax, const vector<double> &u, vector<vector<double>> &points) {

	double tk = 0;
	double tkPlus1 = 0;

	vector<double> yk = u;
	vector<double> ykMinus1 = u;
	vector<double> ykPlus1 = u;

	double tauK, tauKminus1;
	tauK = tauKminus1 = tauMin;

	int N = yk.size();
	for (int i = 0; i < N; i++) {
		cout << setw(13) << "y" << i + 1;
	}
	cout << setw(13) << "tk\n";



	while (tk < T) {
		tkPlus1 = tk + tauK;
		newton(f, yk, ykPlus1, tauK, tkPlus1, 1e-5, 1e-5, 50);

		vector<double> epsK;
		epsK.reserve(f.size());
		for (int i = 0; i < f.size(); i++) {
			epsK.push_back(-tauK / (tauK + tauKminus1) * (ykPlus1[i] - yk[i] - tauK / tauKminus1 * (yk[i] - ykMinus1[i])));
		}

		bool flag = false;
		for (int i = 0; i < epsK.size(); i++) {
			if (fabs(epsK[i]) > epsI) {
				flag = true;
				break;
			}
		}

		if (flag) {
			tauK /= 2;
			tkPlus1 = tk;
			ykPlus1 = yk;
		}
		else {
			double tauKplus1 = calcTauKplus1Three(epsI, epsK, tauK);
			if (tauKplus1 > tauMax) tauKplus1 = tauMax;
			//if (tauKplus1 < tauMin) tauKplus1 = tauMin;

			vector<double> point;
			point.push_back(tkPlus1);
			for (int i = 0; i < N; i++) {
				point.push_back(ykPlus1[i]);
				cout << setw(14) << ykPlus1[i];
			}
			cout << setw(14) << tkPlus1 << endl;

			points.push_back(point);
			ykMinus1 = yk;
			yk = ykPlus1;
			tauKminus1 = tauK;
			tauK = tauKplus1;
			tk = tkPlus1;
		}
	}

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

double function(const vector<func> &f, vector<double> yk, vector<double> ykPlus1,
	double tauK, double tKplus1, int i) {
	return ykPlus1[i] - yk[i] - tauK * f[i](ykPlus1, tKplus1);
}

void calcJInc(vector<vector<double>>& J, const vector<func> &f, vector<double> yk, vector<double> ykPlus1,
	double tauK, double tKplus1) {
	J.clear();
	J.reserve(f.size());
	vector<double> tmp;
	for (int i = 0; i < f.size(); i++) J.push_back(tmp);
	const double M = 0.01;
	for (int i = 0; i < f.size(); i++) {
		double xM = ykPlus1[i] * M;
		for (int j = 0; j < f.size(); j++) {
			double fx = function(f, yk, ykPlus1, tauK, tKplus1, j);
			ykPlus1[i] += xM;
			J[j].push_back((function(f, yk, ykPlus1, tauK, tKplus1, j) - fx) / xM);
			ykPlus1[i] -= xM;
		}
	}
}

int newton(const vector<func> &f, const vector<double>& yk, vector<double> &ykPlus1,
	double tauK, double tKplus1, double e1, double e2, int NIT) {
	int k = 1;
	vector<double> F = yk;
	if (ykPlus1.size() != f.size()) {
		ykPlus1.clear();
		ykPlus1.reserve(f.size());
		for (int i = 0; i < f.size(); i++) {
			ykPlus1.push_back(0.0);
		}
	}
	vector<vector<double>> J;
	vector<double> dx;
	vector<double> resk = ykPlus1;
	double d1, d2;
	do {
		for (int i = 0; i < f.size(); i++) {
			F[i] = -function(f, yk, ykPlus1, tauK, tKplus1, i);
		}
		calcJInc(J, f, yk, ykPlus1, tauK, tKplus1);

		if (gaus(J, F, dx) != 0) {
			cout << "gaus gg\n";
		}

		for (int i = 0; i < resk.size(); i++) {
			resk[i] = ykPlus1[i] + dx[i];
		}

		d1 = DBL_MIN;
		for (int i = 0; i < ykPlus1.size(); i++) {
			double tmp = fabs(function(f, yk, ykPlus1, tauK, tKplus1, i));
			if (tmp > d1) d1 = tmp;
		}

		d2 = DBL_MIN;
		for (int i = 0; i < resk.size(); i++) {
			double tmp = fabs(resk[i] - ykPlus1[i]);
			if (resk[i] >= 1) {
				tmp /= resk[i];
			}

			if (tmp > d2) d2 = tmp;
		}

		ykPlus1 = resk;

		if (k >= NIT) {
			cout << "newton gg " << tKplus1 << endl;
			return 2;
		}
		k++;
	} while (d1 > e1 && d2 > e2);
	//cout << setw(100) << k << endl;
	return 0;
}