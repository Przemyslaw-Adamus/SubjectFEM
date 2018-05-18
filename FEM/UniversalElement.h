#pragma once
#include <cmath>
#include <iostream>

class UniversalElement
{
public:
	double Ksi[4]; //  { (-1 / sqrt(3)),(1 / sqrt(3)),(1 / sqrt(3)),(-1 / sqrt(3)) };
	double Eta[4]; //  { (-1 / sqrt(3)),(-1 / sqrt(3)),(1 / sqrt(3)),(1 / sqrt(3)) };

	double N[4][4];
	double dNdE[4][4];
	double dNdn[4][4];

	double N1N1T[4][4];
	double N2N2T[4][4];
	double N3N3T[4][4];
	double N4N4T[4][4];

	UniversalElement();

	// N1 
	double dN1poKsi(double n);
	double dN1poEta(double n);
	// N2 
	double dN2poKsi(double n);
	double dN2poEta(double n);
	// N3 
	double dN3poKsi(double n);
	double dN3poEta(double n);
	// N4 
	double dN4poKsi(double n);
	double dN4poEta(double n);

	//Funkcje kszta³tu elementu czterowêz³owego
	double N0(double a, double b);
	double N1(double a, double b);
	double N2(double a, double b);
	double N3(double a, double b);

	void print();

	double dxpoKsi(double n, double x1, double x2, double x3, double x4);
	double dxpoEta(double n, double x1, double x2, double x3, double x4);
	double dypoKsi(double n, double y1, double y2, double y3, double y4);
	double dypoEta(double n, double y1, double y2, double y3, double y4);
};

	
