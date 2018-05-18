#pragma once
class MatrixJacobian
{
public:
	double J[2][2];
	double JT[2][2];
	double detJ;
	
	MatrixJacobian();
	double determinant(double J[2][2]);
	void inverse();
};