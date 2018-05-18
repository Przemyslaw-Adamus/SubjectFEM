#pragma once
#include "Node.h"
#include "MatrixJacobian.h"

class Element
{
public:

	int elementId;
	int nodeId[4];
	MatrixJacobian matrixJacobian[4];
	double matrixH[4][4];
	double matrixC[4][4];

	double dNdx[4][4];
	double dNdy[4][4];

	double dNdxdNdxTDetJ1[4][4];
	double dNdxdNdxTDetJ2[4][4];
	double dNdxdNdxTDetJ3[4][4];
	double dNdxdNdxTDetJ4[4][4];

	double dNdydNdyTDetJ1[4][4];
	double dNdydNdyTDetJ2[4][4];
	double dNdydNdyTDetJ3[4][4];
	double dNdydNdyTDetJ4[4][4];

	double N1N1TcRo[4][4];
	double N2N2TcRo[4][4];
	double N3N3TcRo[4][4];
	double N4N4TcRo[4][4];

	double sideLengths[4];
	double boundaryConditionsMatrix[4][4];
	bool boundaryConditions[4];
	double sum1[4][4];
	double sum2[4][4];
	double sum3[4][4];
	double sum4[4][4];

	double dNdxDetJ1[4];
	double dNdxDetJ2[4];
	double dNdxDetJ3[4];
	double dNdxDetJ4[4];

	double dNdyDetJ1[4];
	double dNdyDetJ2[4];
	double dNdyDetJ3[4];
	double dNdyDetJ4[4];

public:
	Element();

};