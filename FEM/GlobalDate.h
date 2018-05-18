#pragma once
#include <cmath>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <iostream>

class GlobalData
{
public:
	std::string fileName;
	std::fstream file;
	double lenghtH;
	double lenghtB;
	int nh;
	int nb;
	double k;
	double c;
	double ro;
	double alfa;
	int dTau;
	double ambientTemperature;
	int simulationTime;
	double initialTemperature;


	GlobalData(); //ustawiamy œcie¿kê do pliku
	bool loadData();// wczytywanie z pliku
	bool saveNetwork();//zapis do pliku
};
