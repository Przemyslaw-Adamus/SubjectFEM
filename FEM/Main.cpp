#include <iostream>
#include "Network.h"

using namespace std;

int main()
{
	Network network;
	UniversalElement universalElement;
	network.createNet();

	std::cout << "___________________________________________" << std::endl;
	std::cout << "   NODES:" << std::endl;
	std::cout << "___________________________________________" << std::endl;
	for (int i = 0;i < network.nodeAmount;i++)
	{
		cout << "ID: "<< network.nodeArray[i].idNode << "\t x: " << network.nodeArray[i].coordX << " \t y "<< network.nodeArray[i].coordY << endl;
	}
	std::cout << "___________________________________________" << std::endl;
	std::cout << "   ELEMENTS:" << std::endl;
	std::cout << "___________________________________________" << std::endl;
	for (int i = 0;i < network.elementsAmount;i++)
	{
		cout << "ID: " << network.elementArray[i].elementId << " ->\t" << network.elementArray[i].nodeId[0] << "\t" << network.elementArray[i].nodeId[1] << "\t" << network.elementArray[i].nodeId[2] << "\t" << network.elementArray[i].nodeId[3] << endl;
	}

	for (int k = 0;k < network.elementsAmount;k++)
	{
		cout << endl << "MATRIX H" << endl;
		for (int i = 0;i < 4;i++)
		{
			for (int j = 0;j < 4;j++)
			{
				cout << network.elementArray[k].matrixH[i][j] << "\t";
			}
			cout << endl;
		}
		cout << endl;
	}

	cout << endl << "MATRIX C" << endl;
	for (int i = 0;i < 4;i++)
	{
		for (int j = 0;j < 4;j++)
		{
			cout << network.elementArray[0].matrixC[i][j] << "\t\t";
		}
		cout << endl;
	}

	for (int k = 0;k < network.elementsAmount;k++)
	{
		cout << endl << "MATRIX Boundary Conditions" << endl;
		for (int i = 0;i < 4;i++)
		{
			for (int j = 0;j < 4;j++)
			{
				cout << network.elementArray[k].boundaryConditionsMatrix[i][j] << "\t\t";
			}
			cout << endl;
		}
	}

	cout << endl << "MATRIX H GENERAL" << endl;
	for (int i = 0;i < network.nodeAmount;i++)
	{
		for (int j = 0;j < network.nodeAmount;j++)
		{
			cout << network.matrixH[i][j] << "\t";
		}
		cout << endl;
	}

	cout << endl << "MATRIX C GENERAL" << endl;
	for (int i = 0;i < network.nodeAmount;i++)
	{
		for (int j = 0;j < network.nodeAmount;j++)
		{
			cout << network.matrixC[i][j] << "\t";
		}
		cout << endl;
	}
	cout << endl;

	network.symulation();

	cout << endl << "MATRIX H+C/dTau GENERAL" << endl;
	for (int i = 0;i < network.nodeAmount;i++)
	{
		for (int j = 0;j < network.nodeAmount;j++)
		{
			cout << network.matrixH[i][j] << "\t";
		}
		cout << endl;
	}
	cout << endl;

	cout << endl << "VECTOR P GENERAL" << endl;
	for (int i = 0;i < network.nodeAmount;i++)
	{
		cout << network.vectorP[i] << "\t";
	}
	system("pause");
	return 0;
}
