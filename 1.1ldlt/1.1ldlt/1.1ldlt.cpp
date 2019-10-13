#include "pch.h"
#include <iomanip>
#include <time.h>
#include <iostream>
using namespace std;

int LDLT(double** A, double* B, int n, double* &res);

void calcF(double* F, double** A, double* B, double* res, int n);
double calcNorm(double* F, int n);
double calcEps(double* res, double* res2, int n);


int main(){
	setlocale(LC_ALL, "rus");
	srand(int(time(NULL)));

	cout << "количество строк и столбцов n = 3\n\n";
	const int n = 3;
	double* B = new double[n];
	double** A = new double*[n];
	for (int i = 0; i < n; i++) {
		A[i] = new double[n];
	}

	cout << "2*l1 + 4*l2    2*(l1 - l2)        2*(l1 - l2) \n";
	cout << "2*(l1 - l2)    2*l1 + l2 + 3*l3   2*l1 + l2 - 3*l3 \n";
	cout << "2*(l1 - l2)    2*l1 + l2 - 3*l3   2*l1 + l2 + 3*l3 \n\n";

	cout << "B = -4*l1 - 2*l2   -4*l1 + l2 + 9*l3   -4*l1 + l2 - 9*l3\n\n";

	cout << "1 - Случайные l1 l2 l3\n2 - Ручной ввод\n";
	int inp;
	while (true) {
		cin >> inp;
		if (inp == 2 || inp == 1) break;
		cout << "Введите 1 или 2\n";
	}
	int l[n];
	if (inp == 2) {
		cout << "l1, l2, l3 = ";
		for (int i = 0; i < n; i++) {
			cin >> l[i];
		}
	}
	else {
		cout << "l1, l2, l3 = ";
		for (int i = 0; i < n; i++) {
			l[i] = rand() % 1000001;
			cout << setw(7) << l[i];
		}
		cout << endl;
	}

	A[0][0] = 2 * l[0] + 4 * l[1];
	A[1][0] = A[2][0] = A[0][1] = A[0][2] = 2 * (l[0] - l[1]);
	A[1][1] = A[2][2] = 2 * l[0] + l[1] + 3 * l[2];
	A[2][1] = A[1][2] = 2 * l[0] + l[1] - 3 * l[2];

	B[0] = -4 * l[0] - 2 * l[1];
	B[1] = -4 * l[0] + l[1] + 9 * l[2];
	B[2] = -4 * l[0] + l[1] - 9 * l[2];


	cout << "A | B\n";
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			cout << setw(14) << A[i][j];
		}
		cout << "  | " << setw(14) << B[i];
		cout << endl;
	}

	double* res = new double[n];
	if (LDLT(A, B, n, res) == 1) {
		cout << "A - не положительно определенная\n";
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
	LDLT(A, B2, n, res2);

	double eps = calcEps(res, res2, n);
	cout << "Eps =   " << setw(14) << scientific << eps << endl;


}



int LDLT(double** A, double* B, int n, double* &res) {

	double** L = new double*[n];
	for (int i = 0; i < n; i++) L[i] = new double[n];
	double* D = new double[n];

	//Разложение матрицы A на L и D
	double sum = 0;
	for (int i = 0; i < n; i++) {
		for (int j = i; j < n; j++) {
			sum = A[j][i];//значение вычисляемого элемента
			for (int k = 0; k < i; k++)//вычитание элементов строки из вычисляемого элемента
				sum = sum - L[i][k] * D[k] * L[j][k];
			if (i == j) {
				if (fabs(sum) <= DBL_EPSILON) {
					for (int i = 0; i < n; i++) delete[] L[i];
					delete[] L;
					delete[] D;
					return 1;//A не положительно определена;
				}
				D[i] = sum;//диагональный элемент
				L[i][i] = 1;//диагональ
			}
			else L[j][i] = sum / D[i];//внедиагональный элемент
		}
	}

	double* Y = new double[n];
	//LY= B  //L - нижнетреугольная матрица с единичн диаг
	for (int i = 0; i < n; i++) {
		Y[i] = B[i];
		for (int k = 0; k < i; k++) {
			Y[i] -= L[i][k] * Y[k];
		}
	}
	// DZ = Y //
	double* Z = new double[n];
	for (int i = 0; i < n; i++) {
		Z[i] = Y[i] / D[i];
	}


	for (int i = n - 1; i >= 0; i--) {
		res[i] = Z[i];
		for (int j = i + 1; j < n; j++) {
			res[i] -= L[j][i] * res[j];
		}
	}

	for (int i = 0; i < n; i++) {
		delete[] L[i];
	}
	delete[] L;
	delete[] D;
	delete[] Y;
	delete[] Z;

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