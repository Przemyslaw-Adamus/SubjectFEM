#include "MatrixJacobian.h"

MatrixJacobian::MatrixJacobian()
{
	J[0][0] = 0.0;
	J[1][0] = 0.0;
	J[0][1] = 0.0;
	J[1][1] = 0.0;

	JT[0][0] = 0.0;
	JT[1][0] = 0.0;
	JT[0][1] = 0.0;
	JT[1][1] = 0.0;

	detJ = 0.0;
}

double MatrixJacobian::determinant(double J[2][2])
{
	return J[0][0] * J[1][1] - J[1][0] * J[0][1];
}

void MatrixJacobian::inverse()
{
	JT[0][0] = J[1][1] / detJ;
	JT[1][1] = J[0][0] / detJ;
	JT[0][1] = -J[0][1] / detJ;
	JT[1][0] = -J[1][0] / detJ;
}