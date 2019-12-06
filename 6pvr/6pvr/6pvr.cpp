#include "pch.h"
#include <iostream>
#include <stdlib.h> 
#include <math.h> 
#include <fstream> 
#define  PRES   double 
#define  NXB	20
#define  NX     NXB*3+1 
#define  NYB    20
#define  NY     NYB*3+1 
#define  REP    30000 
#define  EPSL   1.e-6 
#define  LL     1.7f 
#define  TEM1   20.0f 
#define  TEM2   10.0f 
#define  HX     0.2f 
#define  HY     0.3f 
using namespace std;

void maxpvr(PRES *t1, PRES *del, PRES *maxdel) {
	PRES d = fabs(*del) / fabs(*t1);
	if (d > *maxdel) *maxdel = d;
}
int main(int argc, char **argv) {
	ofstream foutT("dT.dat", ios_base::out | ios_base::trunc | ios_base::binary);
	int   i1, i2, i3, j1, j2, j3, rp, i, j, k = 0;
	PRES  T1 = TEM1, T2 = TEM2, h = HX, r = HY, tx, t0, t1, del, maxdel = 0.0f;
	//PRES  T[NY][NX];
	PRES** T = new PRES*[NY];
	for (int i = 0; i < NY; i++) T[i] = new PRES[NX];
	PRES  lam = LL;
	PRES  eps = EPSL;
	int    prz = 1;

	int    nT = 0;
	PRES  alf_1 = -h / r;
	PRES  alf_2 = -r / h;
	PRES  alf_3 = alf_2 * 0.5f;
	PRES  alf_4 = alf_1 * 0.5f;
	PRES  gam_1 = -2.f * (alf_1 + alf_2);
	PRES  gam_2 = -1.5f * (alf_1 + alf_2);
	PRES  gam_3 = -(alf_1 + alf_2);
	PRES  gam_4 = -(alf_3 + alf_4);
	i1 = NXB; i2 = 2*NXB;   i3 = 3*NXB;
	j1 = NYB; j2 = NYB * 2; j3 = NYB * 3;   rp = REP;
	for (j = 0; j <= j3; j++) {
		for (i = 0; i <= i3; i++) { T[j][i] = 0.0f; }
	}
	for (j = 0; j <= NYB * 3; j++) T[j][0] = T1;
	for (i = 0; i <= i1; i++) T[3*NYB][i] = T1;
	for (j = 2*NYB; j <= 3 * NYB; j++) T[j][3*NXB] = T2;
	while (k < rp && prz == 1) {
		k++;
		for (j = 0; j <= j3; j++) {
			for (i = 1; i <= i3; i++) {
				t0 = T[j][i];
				//2 CD
				if (i > i1 && i < i3 && j == j3) {
					tx = -(alf_3*(T[j][i - 1] + T[j][i + 1]) + alf_1 * T[j - 1][i]) / gam_3;
					del = lam * (tx - t0); t1 = t0 + del; T[j][i] = t1;
					maxpvr(&t1, &del, &maxdel);
				}
				//5.7
				else if (i == i2 && j > j1 && j < j2 || i == i1 && j > 0 && j < j1) {
					tx = -(alf_4*(T[j - 1][i] + T[j + 1][i]) + alf_2*T[j][i - 1])/gam_3;
					del = lam * (tx - t0); t1 = t0 + del; T[j][i] = t1;
					maxpvr(&t1, &del, &maxdel);
				}
				//4.6.8
				else if (i > i2 && i < i3 && j == j2 || i > i1 && i < i2 && j == j1
					|| i > 0 && i < i1 && j == 0) {
					tx = -(alf_3 * (T[j][i - 1] + T[j][i + 1]) + alf_1 * T[j + 1][i]) / gam_3;
					del = lam * (tx - t0); t1 = t0 + del; T[j][i] = t1;
					maxpvr(&t1, &del, &maxdel);
				}
				//10.12.14
				else if (i == i3 && j == j2 || i == i2 && j == j1 || i == i1 && j == 0) {
					tx = -(alf_3*T[j][i - 1] + alf_4 * T[j + 1][i]) / gam_4;
					del = lam * (tx - t0); t1 = t0 + del; T[j][i] = t1;
					maxpvr(&t1, &del, &maxdel);
				}
				//11.13
				
				else if (i == i2 && j == j2 || i == i1 && j == j1) {
					tx = -(alf_4*T[j - 1][i] + alf_2 * T[j][i - 1] + alf_3 * T[j][i + 1] + alf_1 * T[j + 1][i]) / gam_2;
					del = lam * (tx - t0); t1 = t0 + del; T[j][i] = t1;
					maxpvr(&t1, &del, &maxdel);
				}
				//15
				else if (i > 0 && i < i2 && j > j2 && j < j3 || i >= i2 && i < i3 && j > j2 && j <= j3
					|| i > 0 && i < i2 && j > j1 && j <= j2 || i > 0 && i < i1 && j > 0 && j <= j1) {
					tx = -(alf_1*(T[j - 1][i] + T[j + 1][i]) + alf_2 * (T[j][i - 1] + T[j][i + 1])) / gam_1;
					del = lam * (tx - t0); t1 = t0 + del; T[j][i] = t1;
					maxpvr(&t1, &del, &maxdel);
				}
			}
		}
		nT++;  PRES w = maxdel;
		foutT.write((char*)&w, sizeof w);
		if (maxdel < eps) prz = 0; maxdel = 0.0f;
	}
	cout << k << endl;
	foutT.close();
	ofstream fouT("nT.dat", ios_base::out | ios_base::trunc | ios_base::binary);
	fouT.write((char*)&nT, sizeof nT);
	fouT.close();
	ofstream fout("Pole.dat", ios_base::out | ios_base::trunc | ios_base::binary);
	for (j = 0; j < NY; j++) {
		for (i = 0; i < NX; i++) {
			PRES w = T[j][i];
			fout.write((char*)&w, sizeof w);
		}
	}
	fout.close();
	int n_x = NX;  int n_y = NY;
	ofstream fou("Param.dat", ios_base::out | ios_base::trunc | ios_base::binary);
	fou.write((char*)&n_x, sizeof n_x);
	fou.write((char*)&n_y, sizeof n_y);
	fou.close();
	for (int i = 0; i < NY; i++) delete[] T[i];
	delete[] T;
	return 0;
}