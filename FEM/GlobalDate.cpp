#include "GlobalDate.h"


GlobalData::GlobalData()
{
	this->fileName = "data2D.txt";
}

bool GlobalData::loadData()
{
	fileName = "data2D.txt";
	file.open(fileName, std::ios::in);
	if (file.good())
	{
		std::cout << "File was opened succesfully\n";
		file >> initialTemperature;
		file >> simulationTime;
		file >> dTau;
		file >> ambientTemperature;
		file >> alfa;
		file >> lenghtH;
		file >> lenghtB;
		file >> nb;
		file >> nh;
		file >> c;
		file >> k;
		file >> ro;
		
		file.close();
		std::cout << "Data is read\n";
		return true;
	}
	else
	{
		std::cout << "Unable to open file\n";
		file.close();
		return false;
	}
}

bool GlobalData::saveNetwork()
{
	fileName = "Net2D.csv";
	file.open(fileName, std::ios::app);
	if (file.good())
	{
		/*for (int i = 0; i < net.nodeAmount;i++)
		{
			file << net.nodeArray[i].coordX << net.nodeArray[i].coordY << net.nodeArray[i].idNode << std::endl;
		}*/
		return true;
	}
	else
	{
		return false;
	}
	
}
	
