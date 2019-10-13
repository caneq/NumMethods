#include "pch.h"
#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;

template<class T>
void print(T** A, T* B, int n);
int gaus(double** A, double* B, int n, double* &res);
int plot(double* f, double* X, double* Y, int n, int m);

int mnk(double* X, double* Y, int n, int m, double* res, double& disp) {
	int m2 = 2 * m;
	double* POWERX = new double[m2];
	for (int i = 0; i < m2; i++) POWERX[i] = 0;
	for (int i = 0; i < n; i++) {
		double powXi = X[i];
		for (int k = 0; k < m2; k++) {
			POWERX[k] += powXi;
			powXi *= X[i];
		}
	}
	//cout << "POWERX = ";
	//for (int i = 0; i < m2; i++) cout << setw(15) << POWERX[i];
	//cout << endl << endl;

	double** SUMX = new double*[m + 1];
	for (int i = 0; i < m + 1; i++) {
		SUMX[i] = new double[m + 1];
	}

	for (int l = 0; l < m + 1; l++) {
		for (int j = 0; j < m + 1; j++) {
			SUMX[l][j] = l + j - 1 >= 0 ? POWERX[l + j - 1] : n;
		}
	}


	double* PRAW = new double[m + 1];
	for (int i = 0; i < m + 1; i++) PRAW[i] = 0;

	for (int i = 0; i < n; i++) PRAW[0] += Y[i];
	for (int i = 0; i < n; i++) {
		double powXi = X[i];
		for (int l = 1; l < m + 1; l++) {
			PRAW[l] += powXi * Y[i];
			powXi *= X[i];
		}
	}
	//cout << "SUMX | PRAW\n";
	//print(SUMX, PRAW, m + 1);
	gaus(SUMX, PRAW, m + 1, res);


	disp = 0;
	for (int i = 0; i < n; i++) {
		double sum = res[0];
		double powX = X[i];
		for (int j = 1; j < m + 1; j++) {
			sum += res[j] * powX;
			powX *= X[i];
		}
		sum = -sum;
		sum += Y[i];
		disp += sum * sum;
	}
	disp /= n - m - 1;

	return 0;
}

int main() {
	int n;
	cout << "n = ";
	cin >> n;
	if (n <= 0) return 0;
	int m;
	cout << "m = ";
	cin >> m;
	if (m <= 0 || m >= n) return 0;
	double* X = new double[n];
	double* Y = new double[n];
	double* res = new double[m + 1];
	double disp;
	cout << "X = ";
	for (int i = 0; i < n; i++) cin >> X[i];
	cout << "\nY = ";
	for (int i = 0; i < n; i++) cin >> Y[i];
	mnk(X, Y, n, m, res, disp);
	for (int i = 0; i < m + 1; i++) cout << setw(15) << res[i];
	cout << "\n\n    disp           = " << disp << endl;
	cout << "standard deviation = " << sqrt(disp) << endl;
	cout << endl << endl;

	if (fabs(res[0]) > DBL_EPSILON) cout << res[0] << (0 != m ? (res[1] > 0 ? "+" : "") : "");
	for (int i = 1; i < m + 1; i++) {
		if (fabs(res[i]) > DBL_EPSILON) cout << res[i] << "*" << "x^" << i << (i != m ? (res[i + 1] > 0 ? "+" : "") : "");
	}
	cout << endl;
	//for (int i = 0; i < n; i++) {
	//	cout << '(' << X[i] << ',' << Y[i] << ')' << endl;
	//}

	plot(res, X, Y, n, m);
}


template<class T>
void print(T** A, T* B, int n) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			cout << setw(10) << A[i][j];
		}
		cout << setw(10) << B[i] << endl;
	}
	cout << endl;
}

int gaus(double** A, double* B, int n, double* &res) {

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

int plot(double* f, double* X, double* Y, int n, int m) {
	ofstream kv("..\\kv.m");
	if (kv.bad())return 1;
	kv << "x1=[";
	for (int i = 0; i < n - 1; i++) kv << X[i] << ',';
	kv << X[n - 1] << "];\n";
	kv << "y1=[";
	for (int i = 0; i < n - 1; i++) kv << Y[i] << ',';
	kv << Y[n - 1] << "];\n";
	kv << "x = (min(x1)-100):0.01:(max(x1)+100);\n";
	kv << "y =";
	if (fabs(f[0]) > DBL_EPSILON) kv << f[0] << (0 != m ? (f[1] > 0 ? "+" : "") : "");
	for (int i = 1; i < m + 1; i++) {
		if (fabs(f[i]) > DBL_EPSILON) kv << f[i] << ".*" << "x.^" << i << (i != m ? (f[i + 1] > 0 ? "+" : "") : "");
	}
	kv << ";\n";
	kv << "plot(x,y,x1,y1,'r.','MarkerSize', 15);\n";
	kv << "xlim([min(x1),max(x1)]);\n";
	kv << "ylim([min(y1),max(y1)]);\n";
	kv << "grid on;\n";
	kv << "xlabel('t, градусы C');\n";
	kv << "ylabel('r, мкОм');\n";
}