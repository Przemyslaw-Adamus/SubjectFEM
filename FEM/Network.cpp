#include "Network.h"
#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>

Network::Network()
{
	loadData();
	nodeAmount = nb * nh;
	elementsAmount = (nh - 1)*(nb - 1);

	elementArray = new Element[elementsAmount]();
	nodeArray = new Node[nodeAmount]();
	matrixH = new double*[nodeAmount];
	matrixC = new double*[nodeAmount];
	matrixHGeneral = new double*[nodeAmount];
	vectorP = new double[nodeAmount];

	for (int i = 0;i < nodeAmount;i++)
	{
		matrixH[i] = new double[nodeAmount];
		matrixC[i] = new double[nodeAmount];
		matrixHGeneral[i] = new double[nodeAmount];
	}

	for (int i = 0;i < nodeAmount;i++)
	{
		for (int j = 0;j < nodeAmount;j++)
		{
			matrixH[i][j] = 0.0;
			matrixC[i][j] = 0.0;
			matrixHGeneral[i][j] = 0.0;
		}
		vectorP[i] = 0.0;
	}

}

void Network::createNet()
{
	/*globalData.loadData();
	nodeAmount = globalData.nb*globalData.nh;
	elementsAmount = (globalData.nh - 1)*(globalData.nb - 1);*/
	Network::setNodes();
	Network::setElemens();
	Network::setJacobianMatrix();
	Network::createHMatrix();
	Network::createCMatrix();
	Network::setSideLengths();
	Network::setBoundaryConditions();
	Network::setBoundaryConditionsMatrix();
	Network::createP();
}

void Network::setNodes() //Ustawiamy wêz³y
{
	int counterNode = 1;
	for (int i = 1;i < nodeAmount;i++)
	{
		nodeArray[i].idNode = i;
		counterNode = i % nb;
		if (counterNode == 0)
		{
			nodeArray[i].coordX = nodeArray[i - 1].coordX + (lenghtH / (nh - 1));
		}
		else
		{
			nodeArray[i].coordX = nodeArray[i - 1].coordX;
		}
		nodeArray[i].coordY = (lenghtB / (nb - 1))*counterNode;
	}
}

void Network::setElemens() //Ustawiamy elementy
{
	
	for (int i = 0;i < nh - 1;i++)
	{
		for (int j = 0;j < nb - 1;j++)
		{
			elementArray[i * (nb - 1) + j].elementId = i * (nb - 1) + j;
			elementArray[i * (nb - 1) + j].nodeId[0] = nodeArray[i + i * (nb - 1) + j].idNode;
			elementArray[i * (nb - 1) + j].nodeId[1] = nodeArray[i + i * (nb - 1) + j + nb].idNode;
			elementArray[i * (nb - 1) + j].nodeId[2] = nodeArray[i + i * (nb - 1) + j + nb + 1].idNode;
			elementArray[i * (nb - 1) + j].nodeId[3] = nodeArray[i + i * (nb - 1) + j + 1].idNode;
		}
	}
}

void Network::setJacobianMatrix() //Ustawiamy macierze Jacobian 4*ne
{
	for (int i = 0; i < elementsAmount;i++)
	{
		for (int j = 0;j < 4;j++)
		{
			elementArray[i].matrixJacobian[j].J[0][0] =
				dxpoKsi(Ksi[j], nodeArray[elementArray[i].nodeId[0]].coordX, nodeArray[elementArray[i].nodeId[1]].coordX, nodeArray[elementArray[i].nodeId[2]].coordX, nodeArray[elementArray[i].nodeId[3]].coordX);
			elementArray[i].matrixJacobian[j].J[0][1] =
				dypoKsi(Ksi[j], nodeArray[elementArray[i].nodeId[0]].coordY, nodeArray[elementArray[i].nodeId[1]].coordY, nodeArray[elementArray[i].nodeId[2]].coordY, nodeArray[elementArray[i].nodeId[3]].coordY);
			elementArray[i].matrixJacobian[j].J[1][0] =
				dxpoEta(Eta[j], nodeArray[elementArray[i].nodeId[0]].coordX, nodeArray[elementArray[i].nodeId[1]].coordX, nodeArray[elementArray[i].nodeId[2]].coordX, nodeArray[elementArray[i].nodeId[3]].coordX);
			elementArray[i].matrixJacobian[j].J[1][1] =
				dypoEta(Eta[j], nodeArray[elementArray[i].nodeId[0]].coordY, nodeArray[elementArray[i].nodeId[1]].coordY, nodeArray[elementArray[i].nodeId[2]].coordY, nodeArray[elementArray[i].nodeId[3]].coordY);

			elementArray[i].matrixJacobian[j].detJ = elementArray[i].matrixJacobian[j].determinant(elementArray[i].matrixJacobian[j].J);
			elementArray[i].matrixJacobian[j].inverse();
		}
	}
}

void Network::createHMatrix()
{
	//Ustawiamy macierze dNdy, dNdx
	for (int i = 0;i < elementsAmount;i++)
	{
		for (int j = 0;j < 4;j++)
		{
			elementArray[i].dNdy[j][0] = elementArray[i].matrixJacobian[j].JT[1][0] * dNdE[0][j] + elementArray[i].matrixJacobian[j].JT[1][1] * dNdn[0][j];
			elementArray[i].dNdy[j][1] = elementArray[i].matrixJacobian[j].JT[1][0] * dNdE[1][j] + elementArray[i].matrixJacobian[j].JT[1][1] * dNdn[1][j];
			elementArray[i].dNdy[j][2] = elementArray[i].matrixJacobian[j].JT[1][0] * dNdE[2][j] + elementArray[i].matrixJacobian[j].JT[1][1] * dNdn[2][j];
			elementArray[i].dNdy[j][3] = elementArray[i].matrixJacobian[j].JT[1][0] * dNdE[3][j] + elementArray[i].matrixJacobian[j].JT[1][1] * dNdn[3][j];

			elementArray[i].dNdx[j][0] = elementArray[i].matrixJacobian[j].JT[0][0] * dNdE[0][j] + elementArray[i].matrixJacobian[j].JT[0][1] * dNdn[0][j];
			elementArray[i].dNdx[j][1] = elementArray[i].matrixJacobian[j].JT[0][0] * dNdE[1][j] + elementArray[i].matrixJacobian[j].JT[0][1] * dNdn[1][j];
			elementArray[i].dNdx[j][2] = elementArray[i].matrixJacobian[j].JT[0][0] * dNdE[2][j] + elementArray[i].matrixJacobian[j].JT[0][1] * dNdn[2][j];
			elementArray[i].dNdx[j][3] = elementArray[i].matrixJacobian[j].JT[0][0] * dNdE[3][j] + elementArray[i].matrixJacobian[j].JT[0][1] * dNdn[3][j];
		}
	}

	//Ustawiamy macierze dNdxdNdxTDetJ...
	for (int i = 0;i < elementsAmount;i++)
	{
		for (int j = 0;j < 4;j++)
		{
			for (int k = 0;k < 4;k++)
			{
				elementArray[i].dNdxdNdxTDetJ1[j][k] = elementArray[i].dNdx[j][0] * elementArray[i].dNdx[k][0] * elementArray[i].matrixJacobian[0].detJ;
				elementArray[i].dNdxdNdxTDetJ2[j][k] = elementArray[i].dNdx[j][1] * elementArray[i].dNdx[k][1] * elementArray[i].matrixJacobian[1].detJ;
				elementArray[i].dNdxdNdxTDetJ3[j][k] = elementArray[i].dNdx[j][2] * elementArray[i].dNdx[k][2] * elementArray[i].matrixJacobian[2].detJ;
				elementArray[i].dNdxdNdxTDetJ4[j][k] = elementArray[i].dNdx[j][3] * elementArray[i].dNdx[k][3] * elementArray[i].matrixJacobian[3].detJ;

				elementArray[i].dNdydNdyTDetJ1[j][k] = elementArray[i].dNdy[j][0] * elementArray[i].dNdy[k][0] * elementArray[i].matrixJacobian[0].detJ;;
				elementArray[i].dNdydNdyTDetJ2[j][k] = elementArray[i].dNdy[j][1] * elementArray[i].dNdy[k][1] * elementArray[i].matrixJacobian[1].detJ;
				elementArray[i].dNdydNdyTDetJ3[j][k] = elementArray[i].dNdy[j][2] * elementArray[i].dNdy[k][2] * elementArray[i].matrixJacobian[2].detJ;
				elementArray[i].dNdydNdyTDetJ4[j][k] = elementArray[i].dNdy[j][3] * elementArray[i].dNdy[k][3] * elementArray[i].matrixJacobian[3].detJ;
			}
		}
	}

	//Ustawiamy macierz H
	for (int i = 0;i < elementsAmount;i++)
	{
		for (int j = 0;j < 4;j++)
		{
			for (int m = 0;m < 4;m++)
			{
				matrixH[elementArray[i].nodeId[j]][elementArray[i].nodeId[m]] +=
					k * (elementArray[i].dNdxdNdxTDetJ1[j][m] + elementArray[i].dNdxdNdxTDetJ2[j][m] + elementArray[i].dNdxdNdxTDetJ3[j][m] + elementArray[i].dNdxdNdxTDetJ4[j][m] + elementArray[i].dNdydNdyTDetJ1[j][m] + elementArray[i].dNdydNdyTDetJ2[j][m] + elementArray[i].dNdydNdyTDetJ3[j][m] + elementArray[i].dNdydNdyTDetJ4[j][m]);
				elementArray[i].matrixH[j][m] =
					k*(elementArray[i].dNdxdNdxTDetJ1[j][m] + elementArray[i].dNdxdNdxTDetJ2[j][m] + elementArray[i].dNdxdNdxTDetJ3[j][m] + elementArray[i].dNdxdNdxTDetJ4[j][m] + elementArray[i].dNdydNdyTDetJ1[j][m] + elementArray[i].dNdydNdyTDetJ2[j][m] + elementArray[i].dNdydNdyTDetJ3[j][m] + elementArray[i].dNdydNdyTDetJ4[j][m]);
			}
		}
	}
}

void Network::createCMatrix()
{
	//Ustawiamy macierz N1N1T,N2N2T,N3N3T,N4N4T
	for (int i = 0;i < elementsAmount;i++)
	{
		for (int j = 0;j < 4;j++)
		{
			for (int k = 0;k < 4;k++)
			{
				elementArray[i].N1N1TcRo[j][k] =
					N1N1T[j][k] * c*ro*elementArray[i].matrixJacobian[0].detJ;
				elementArray[i].N2N2TcRo[j][k] =
					N2N2T[j][k] * c*ro*elementArray[i].matrixJacobian[1].detJ;
				elementArray[i].N3N3TcRo[j][k] =
					N3N3T[j][k] * c*ro*elementArray[i].matrixJacobian[2].detJ;
				elementArray[i].N4N4TcRo[j][k] =
					N4N4T[j][k] * c*ro*elementArray[i].matrixJacobian[3].detJ;
			}
		}
	}

	//Ustawiamy macierz C
	for (int i = 0;i < elementsAmount;i++)
	{
		for (int j = 0;j < 4;j++)
		{
			for (int k = 0;k < 4;k++)
			{
				elementArray[i].matrixC[j][k] =
					elementArray[i].N1N1TcRo[j][k] + elementArray[i].N2N2TcRo[j][k] + elementArray[i].N3N3TcRo[j][k] + elementArray[i].N4N4TcRo[j][k];
				matrixC[elementArray[i].nodeId[j]][elementArray[i].nodeId[k]] +=
					elementArray[i].N1N1TcRo[j][k] + elementArray[i].N2N2TcRo[j][k] + elementArray[i].N3N3TcRo[j][k] + elementArray[i].N4N4TcRo[j][k];


			}
		}
	}
}

void Network::setSideLengths()
{
	for (int k = 0;k < elementsAmount; k++)
	{
		for (int i = 0;i < 4;i++)
		{
			elementArray[k].sideLengths[0] = sqrt(pow((nodeArray[elementArray[k].nodeId[0]].coordX - nodeArray[elementArray[k].nodeId[1]].coordX), 2) + pow((nodeArray[elementArray[k].nodeId[0]].coordY - nodeArray[elementArray[k].nodeId[1]].coordY), 2));
			elementArray[k].sideLengths[1] = sqrt(pow((nodeArray[elementArray[k].nodeId[1]].coordX - nodeArray[elementArray[k].nodeId[2]].coordX), 2) + pow((nodeArray[elementArray[k].nodeId[1]].coordY - nodeArray[elementArray[k].nodeId[2]].coordY), 2));
			elementArray[k].sideLengths[2] = sqrt(pow((nodeArray[elementArray[k].nodeId[2]].coordX - nodeArray[elementArray[k].nodeId[3]].coordX), 2) + pow((nodeArray[elementArray[k].nodeId[2]].coordY - nodeArray[elementArray[k].nodeId[3]].coordY), 2));
			elementArray[k].sideLengths[3] = sqrt(pow((nodeArray[elementArray[k].nodeId[3]].coordX - nodeArray[elementArray[k].nodeId[0]].coordX), 2) + pow((nodeArray[elementArray[k].nodeId[3]].coordY - nodeArray[elementArray[k].nodeId[0]].coordY), 2));
		}
		//std::cout << elementArray[k].sideLengths[0] << "\t" << elementArray[k].sideLengths[1] << "\t" << elementArray[k].sideLengths[2] << "\t" << elementArray[k].sideLengths[3] << std::endl;
	}	
}

void Network::setBoundaryConditionsMatrix()
{
	double pc1[4];
	double pc2[4];

	pc1[0] = N0(-1 / sqrt(3), -1);
	pc1[1] = N1(-1 / sqrt(3), -1);
	pc1[2] = N2(-1 / sqrt(3), -1);
	pc1[3] = N3(-1 / sqrt(3), -1);

	pc2[0] = N0(1 / sqrt(3), -1);
	pc2[1] = N1(1 / sqrt(3), -1);
	pc2[2] = N2(1 / sqrt(3), -1);
	pc2[3] = N3(1 / sqrt(3), -1);

	for (int m = 0;m < elementsAmount;m++)
	{
		for (int i = 0;i < 4;i++)
		{
			for (int j = 0;j < 4;j++)
			{
				elementArray[m].sum1[i][j] = ((pc1[i] * pc1[j]) + (pc2[i] * pc2[j]))*elementArray[m].sideLengths[0] / 2;
			}
		}
	}

	pc1[0] = N0(1,-1 / sqrt(3));
	pc1[1] = N1(1,-1 / sqrt(3));
	pc1[2] = N2(1,-1 / sqrt(3));
	pc1[3] = N3(1,-1 / sqrt(3));

	pc2[0] = N0(1,1 / sqrt(3));
	pc2[1] = N1(1,1 / sqrt(3));
	pc2[2] = N2(1,1 / sqrt(3));
	pc2[3] = N3(1,1 / sqrt(3));

	for (int m = 0;m < elementsAmount;m++)
	{
		for (int i = 0;i < 4;i++)
		{
			for (int j = 0;j < 4;j++)
			{
				elementArray[m].sum2[i][j] = ((pc1[i] * pc1[j]) + (pc2[i] * pc2[j]))*elementArray[m].sideLengths[1] / 2;
			}
		}
	}

	pc1[0] = N0(1 / sqrt(3),1);
	pc1[1] = N1(1 / sqrt(3),1);
	pc1[2] = N2(1 / sqrt(3),1);
	pc1[3] = N3(1 / sqrt(3),1);

	pc2[0] = N0(-1 / sqrt(3),1);
	pc2[1] = N1(-1 / sqrt(3),1);
	pc2[2] = N2(-1 / sqrt(3),1);
	pc2[3] = N3(-1 / sqrt(3),1);

	for (int m = 0;m < elementsAmount;m++)
	{
		for (int i = 0;i < 4;i++)
		{
			for (int j = 0;j < 4;j++)
			{
				elementArray[m].sum3[i][j] = ((pc1[i] * pc1[j]) + (pc2[i] * pc2[j]))*elementArray[m].sideLengths[2] / 2;
			}
		}
	}

	pc1[0] = N0(-1, 1 / sqrt(3));
	pc1[1] = N1(-1, 1 / sqrt(3));
	pc1[2] = N2(-1, 1 / sqrt(3));
	pc1[3] = N3(-1, 1 / sqrt(3));

	pc2[0] = N0(-1, -1 / sqrt(3));
	pc2[1] = N1(-1, -1 / sqrt(3));
	pc2[2] = N2(-1, -1 / sqrt(3));
	pc2[3] = N3(-1, -1 / sqrt(3));

	for (int m = 0;m < elementsAmount;m++)
	{
		for (int i = 0;i < 4;i++)
		{
			for (int j = 0;j < 4;j++)
			{
				elementArray[m].sum4[i][j] = ((pc1[i] * pc1[j]) + (pc2[i] * pc2[j]))*elementArray[m].sideLengths[3] / 2;
			}
		}
	}

	for (int m = 0;m < elementsAmount;m++)
	{
		for (int i = 0;i < 4;i++)
		{
			for (int j = 0;j < 4;j++)
			{
				elementArray[m].boundaryConditionsMatrix[i][j] = elementArray[m].sum1[i][j] * elementArray[m].boundaryConditions[0] + elementArray[m].sum2[i][j] * elementArray[m].boundaryConditions[1] + elementArray[m].sum3[i][j] * elementArray[m].boundaryConditions[2] + elementArray[m].sum4[i][j] * elementArray[m].boundaryConditions[3];
				elementArray[m].boundaryConditionsMatrix[i][j] = elementArray[m].boundaryConditionsMatrix[i][j] * alfa;
			}
		}
	}
}

void Network::setBoundaryConditions()
{
	elementArray[0].boundaryConditions[0] = true;
	elementArray[0].boundaryConditions[3] = true;
	elementArray[1].boundaryConditions[3] = true;
	elementArray[2].boundaryConditions[3] = true;
	elementArray[2].boundaryConditions[2] = true;
	elementArray[3].boundaryConditions[0] = true;
	elementArray[5].boundaryConditions[2] = true;
	elementArray[6].boundaryConditions[0] = true;
	elementArray[6].boundaryConditions[1] = true;
	elementArray[7].boundaryConditions[1] = true;
	elementArray[8].boundaryConditions[1] = true;
	elementArray[8].boundaryConditions[2] = true;
}

bool Network::saveNetwork()
{
	fileName = "Net2D.csv";
	file.open(fileName, std::ios::out);
	if (file.good())
	{
		for (int i = 0; i < nodeAmount;i++)
		{
			file << nodeArray[i].coordX <<"\n"<< nodeArray[i].coordY << "\n" << nodeArray[i].idNode << ";" << std::endl;
		}
		return true;
	}
	else
	{
		return false;
	}

}

void Network::symulation()
{
	Network::creataHMatrixGeneral();
	double time = 0.0;
	double minTemp = nodeArray[0].temperature;
	double maxTemp = nodeArray[0].temperature;
	while (time <= simulationTime)
	{
		for (int i = 0; i < nodeAmount;i++)
		{
			//vectorP[i] *= nodeArray[i].temperature;
			nodeArray[i].temperature = nodeArray[i].temperature;
			//....
			if (minTemp > nodeArray[i].temperature)
			{
				minTemp = nodeArray[i].temperature;
			}
			if (maxTemp < nodeArray[i].temperature)
			{
				maxTemp = nodeArray[i].temperature;
			}
			std::cout << nodeArray[i].temperature << "\t";
		}
		std::cout << std::endl << "__________________________________________" << std::endl;
		std::cout << "TIME: " << time << "     MIN: " << minTemp << "    MAX: " << maxTemp << std::endl;
		std::cout << "__________________________________________" << std::endl;
		time += dTau;
	}
}

void Network::creataHMatrixGeneral()
{
	for (int i = 0;i < nodeAmount;i++)
	{
		for (int j = 0;j < nodeAmount;j++)
		{
			matrixH[i][j] += matrixC[i][j] / dTau;
		}
	}

	for (int i = 0;i < elementsAmount;i++)
	{
		for (int j = 0;j < 4;j++)
		{
			for (int m = 0;m < 4;m++)
			{
				matrixH[elementArray[i].nodeId[j]][elementArray[i].nodeId[m]] += elementArray[i].boundaryConditionsMatrix[j][m];
			}
		}
	}
}

void Network::createP()
{
	for (int i = 0;i < elementsAmount;i++)
	{
		for (int j = 0;j < 4;j++)
		{

				elementArray[i].dNdxDetJ1[j] = elementArray[i].dNdx[0][0] * elementArray[i].matrixJacobian[0].detJ*alfa*ambientTemperature;
				elementArray[i].dNdxDetJ2[j] = elementArray[i].dNdx[0][1] * elementArray[i].matrixJacobian[1].detJ*alfa*ambientTemperature;
				elementArray[i].dNdxDetJ3[j] = elementArray[i].dNdx[0][2] * elementArray[i].matrixJacobian[2].detJ*alfa*ambientTemperature;
				elementArray[i].dNdxDetJ4[j] = elementArray[i].dNdx[0][3] * elementArray[i].matrixJacobian[3].detJ*alfa*ambientTemperature;

				elementArray[i].dNdyDetJ1[j] = elementArray[i].dNdy[0][0] * elementArray[i].matrixJacobian[0].detJ*alfa*ambientTemperature;
				elementArray[i].dNdyDetJ2[j] = elementArray[i].dNdy[0][1] * elementArray[i].matrixJacobian[1].detJ*alfa*ambientTemperature;
				elementArray[i].dNdyDetJ3[j] = elementArray[i].dNdy[0][2] * elementArray[i].matrixJacobian[2].detJ*alfa*ambientTemperature;
				elementArray[i].dNdyDetJ4[j] = elementArray[i].dNdy[0][3] * elementArray[i].matrixJacobian[3].detJ*alfa*ambientTemperature;
		}
	}

	double tmp[16] = { 0.0 };
	for (int i = 0;i < elementsAmount;i++)
	{
		for (int j = 0;j < 4;j++)
		{
			//vectorP[elementArray[i].nodeId[j]] +=
			tmp[elementArray[i].nodeId[j]] =(elementArray[i].dNdxDetJ1[j] + elementArray[i].dNdxDetJ2[j] + elementArray[i].dNdxDetJ3[j] + elementArray[i].dNdxDetJ4[j] + elementArray[i].dNdyDetJ1[j] + elementArray[i].dNdyDetJ2[j] + elementArray[i].dNdyDetJ3[j] + elementArray[i].dNdyDetJ4[j]);
		}
	}
	for (int i = 0;i < nodeAmount;i++)
	{
			std::cout << tmp[i] << std::endl;
	}
	for (int i = 0;i < elementsAmount;i++)
	{
		vectorP[i] = -vectorP[i];
	}

	for (int i = 0;i < nodeAmount;i++)
	{
		for (int j = 0;j < nodeAmount;j++)
		{
			vectorP[i] += matrixC[i][j] / dTau;
		}
		vectorP[i] *= 100;
	}
}

