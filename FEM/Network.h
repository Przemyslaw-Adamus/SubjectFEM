#pragma once
#include "Element.h"
#include "Node.h"
#include "GlobalDate.h"
#include "UniversalElement.h"


class Network : public UniversalElement, public GlobalData
{
public:
	int nodeAmount;
	int elementsAmount;
	Element* elementArray;
	Node* nodeArray;
	double** matrixH;
	double** matrixHGeneral;
	double** matrixC;
	double* vectorP;

public:
	Network();
	void createNet();
	void setNodes();
	void setElemens();
	void setJacobianMatrix();
	void createHMatrix();
	void createCMatrix();
	void setSideLengths();
	void setBoundaryConditionsMatrix();
	void setBoundaryConditions();
	bool saveNetwork();
	void symulation();
	void creataHMatrixGeneral();
	void createP();
};