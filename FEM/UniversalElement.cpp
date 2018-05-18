#include "UniversalElement.h"
using namespace std;

UniversalElement::UniversalElement()
{
	Ksi[0] = -1 / sqrt(3);
	Ksi[1] = 1 / sqrt(3);
	Ksi[2] = 1 / sqrt(3);
	Ksi[3] = -1 / sqrt(3);

	Eta[0] = -1 / sqrt(3);
	Eta[1] = -1 / sqrt(3);
	Eta[2] = 1 / sqrt(3);
	Eta[3] = 1 / sqrt(3);

	for (int i = 0; i<4; i++)
	{
			N[i][0] = N0(Ksi[i], Eta[i]);
			N[i][1] = N1(Ksi[i], Eta[i]);
			N[i][2] = N2(Ksi[i], Eta[i]);
			N[i][3] = N3(Ksi[i], Eta[i]);
	}

	for (int i = 0; i<4; i++)
	{
		dNdE[i][0] = dN1poKsi(Eta[i]);
		dNdn[i][0] = dN1poEta(Ksi[i]);
		dNdE[i][1] = dN2poKsi(Eta[i]);
		dNdn[i][1] = dN2poEta(Ksi[i]);
		dNdE[i][2] = dN3poKsi(Eta[i]);
		dNdn[i][2] = dN3poEta(Ksi[i]);
		dNdE[i][3] = dN4poKsi(Eta[i]);
		dNdn[i][3] = dN4poEta(Ksi[i]);
	}

	for (int i = 0;i < 4;i++)
	{
		for (int j = 0;j < 4;j++)
		{
			N1N1T[i][j] = N[i][0] * N[j][0];
			N2N2T[i][j] = N[i][1] * N[j][1];
			N3N3T[i][j] = N[i][2] * N[j][2];
			N4N4T[i][j] = N[i][3] * N[j][3];
		}
	}
}

// N1 
double UniversalElement:: dN1poKsi(double n)
{
	return ((-1.0 / 4.0)*(1.0 - n));
}
double UniversalElement:: dN1poEta(double n)
{
	return ((-1.0 / 4.0)*(1.0 - n));
}
// N2 
double UniversalElement:: dN2poKsi(double n)
{
	return ((1.0 / 4.0)*(1.0 - n));
}
double UniversalElement:: dN2poEta(double n)
{
	return ((-1.0 / 4.0)*(1.0 + n));
}
// N3 
double UniversalElement:: dN3poKsi(double n)
{
	return ((1.0 / 4.0)*(1.0 + n));
}
double UniversalElement:: dN3poEta(double n)
{
	return ((1.0 / 4.0)*(1.0 + n));
}
// N4 
double UniversalElement:: dN4poKsi(double n)
{
	return ((-1.0 / 4.0)*(1.0 + n));
}
double UniversalElement:: dN4poEta(double n)
{
	return ((1.0 / 4.0)*(1.0 - n));
}

//Funkcje kszta³tu elementu czterowêz³owego
double UniversalElement :: N0(double a, double b)
{
	return ((1.0 / 4.0)*(1.0 - a)*(1.0 - b));
}
double UniversalElement :: N1(double a, double b)
{
	return ((1.0 / 4.0)*(1.0 + a)*(1.0 - b));
}
double UniversalElement :: N2(double a, double b)
{
	return ((1.0 / 4.0)*(1.0 + a)*(1.0 + b));
}
double UniversalElement :: N3(double a, double b)
{
	return ((1.0 / 4.0)*(1.0 - a)*(1.0 + b));
}

void UniversalElement :: print() {
	cout << "\n MATRIX N \n";
	for (int i = 0; i<4; i++) {
		for (int j = 0; j<4; j++) {
			cout << N[i][j]<<"\t";
		}
		cout << "\n";
	}
	cout << "\n MATRIX dN/dn \n";
	for (int i = 0; i<4; i++) {
		for (int j = 0; j<4; j++) {
			cout << (dNdn[i][j]) << "\t";
		}
		cout << "\n";
	}
	cout << "\n MATRIX dN/dE \n";
	for (int i = 0; i<4; i++) {
		for (int j = 0; j<4; j++) {
			cout << (dNdE[i][j]) << "\t";
		}
		cout << "\n";
	}
	cout << "\n MATRIX N*NT\n";
	for (int i = 0; i<4; i++) {
		for (int j = 0; j<4; j++) {
			cout << (N1N1T[i][j]) << "\t";
		}
		cout << "\n";
	}
}

double UniversalElement::dxpoKsi(double n, double x1, double x2, double x3, double x4)
{
	return 0.25*(n - 1)*x1 + 0.25*(1 - n)*x2 + 0.25*(1 + n)*x3 - 0.25*(1 + n)*x4;
}
double UniversalElement::dxpoEta(double n, double x1, double x2, double x3, double x4)
{
	return 0.25*(n - 1)*x1 - 0.25*(1 + n)*x2 + 0.25*(1 + n)*x3 + 0.25*(1 - n)*x4;
}
double UniversalElement::dypoKsi(double n, double y1, double y2, double y3, double y4)
{
	return 0.25*(n - 1)*y1 + 0.25*(1 - n)*y2 + 0.25*(1 + n)*y3 - 0.25*(1 + n)*y4;
}
double UniversalElement::dypoEta(double n, double y1, double y2, double y3, double y4)
{
	return 0.25*(n - 1)*y1 - 0.25*(1 + n)*y2 + 0.25*(1 + n)*y3 + 0.25*(1 - n)*y4;
}

