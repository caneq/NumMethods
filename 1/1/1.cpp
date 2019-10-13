#include "pch.h"
#include <iostream>
#include <time.h>
#include <iomanip>
using namespace std;

template<class T>
void print(T** A, T* B, int n);
int gaus(double** A, double* B, int n, double* &res);
void calcF(double* F, double** A, double* B, double* res, int n);
double calcNorm(double* F, int n);
double calcEps(double* res, double* res2, int n);

int main(){
	setlocale(LC_ALL, "rus");
	srand(unsigned int(time(NULL)));

	cout << "количество строк и столбцов n = "; 

	int n;
	while (true) {
		cin >> n;
		if (n > 0) break;
		cout << "n > 0\n";
	}

	double* B = new double[n];
	double** A = new double*[n];
	for (int i = 0; i < n; i++) {
		A[i] = new double[n];
	}

	
	cout << "1 - Случайные числа\n2 - Ручной ввод\n";
	int inp;
	while (true) {
		cin >> inp;
		if (inp == 2 || inp == 1) break;
		cout << "Введите 1 или 2\n";
	}

	if (inp == 2) {
		cout << "A = \n";
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				cin >> A[i][j];
			}
		}

		cout << "B =     ";
		for (int i = 0; i < n; i++) {
			cin >> B[i];
		}
	}
	else {
		cout << "A = \n";
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				A[i][j] = double(rand() % 20001) / 1e2 - 100;
				cout << setw(15) << A[i][j];
			}
			cout << endl;
		}

		cout << "\nB =     ";
		for (int i = 0; i < n; i++) {
			B[i] = double(rand() % 20001) / 1e2 - 100;
			cout << setw(14) << left << B[i];
		}
		cout << endl;
	}
	double* res = new double[n];
	if (gaus(A, B, n, res) == 1) {
		cout << "Решений нет или их бесконечно много\n";
		delete[] A;
		delete[] B;
		delete[] res;
		return 0;
	}
	else {
		cout << "x =     ";
		for (int i = 0; i < n; i++) {
			cout << setw(14) << left << res[i];
		}
		cout << endl;
	}

	double* F = new double[n];
	calcF(F, A, B, res, n);
	cout << "F =     "; 
	for (int i = 0; i < n; i++) cout << setw(14) << scientific << F[i]; cout << endl;

	double norm = calcNorm(F, n);
	cout << "norma = " << setw(14) << scientific << norm << endl;

	double* B2 = new double[n];
	for (int i = 0; i < n; i++) {
		B2[i] = 0;
		for (int j = 0; j < n; j++) {
			B2[i] += res[j] * A[i][j];
		}
	}

	double* res2 = new double[n];
	gaus(A, B2, n, res2);

	double eps = calcEps(res, res2, n);
	cout << "Eps =   " << setw(14) << scientific << eps << endl;

	for (int i = 0; i < n; i++)	delete[] A[i];
	delete[] A;
	delete[] B;
	delete[] res;
	delete[] res2;
	delete[] F;
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

void calcF(double* F, double** A, double* B, double* res, int n) {
	for (int i = 0; i < n; i++) {
		F[i] = 0;
		for (int j = 0; j < n; j++) {
			F[i] += res[j] * A[i][j];
		}
		F[i] -= B[i];
	}
}

double calcNorm(double* F, int n) {
	double norm = fabs(F[0]);
	for (int i = 1; i < n; i++) {
		if (fabs(F[i]) > norm) norm = fabs(F[i]);
	}
	return norm;
}

double calcEps(double* res, double* res2, int n) {
	double eps = fabs(res2[0] - res[0]);
	for (int i = 1; i < n; i++) {
		double tmp = fabs(res2[i] - res[i]);
		if (tmp > eps) eps = tmp;
	}
	double tmp = fabs(res[0]);
	for (int i = 1; i < n; i++) {
		if (fabs(res[i]) > tmp) tmp = fabs(res[i]);
	}
	eps /= tmp;
	return eps;
}