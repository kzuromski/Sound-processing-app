#pragma once
#include "windows.h"
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>
using namespace std;

class Wave
{
private:
	FOURCC riff; // "RIFF" file descripton header
	DWORD size_of_file; // size of file
	FOURCC wave; // "WAVE" file descripton header
	FOURCC fmt; // "fmt" description header
	DWORD chunk; // size of WAVE section chunck
	WORD pcm; // WAVE type format
	WORD chanel; // mono/stereo
	DWORD sample_rate; // sample rate
	DWORD bytes_per_sec; // bytes/sec
	WORD block_alignment; // Block alignment
	WORD bits_per_sample; // Bits/sample
	FOURCC data; // "data" description header
	DWORD size_of_data; // size of data chunk
	const double eps = 1e-12; //wspolczynnik
	INT16 data_of_file;
	vector<INT16> left;
	vector<double> left_minus;
	vector<INT16> right;
	vector<double> right_minus;
	vector<double>minLsrVector;
	ofstream new_file;
	FILE *wf;
	double ammount_of_samples;
	double averageEPS;
	double averageBit;
	double averageLsr;
	int r;

public:
	Wave(string name_of_wave, int r);
	~Wave();
	double getAverageEPS();
	double getAverageBit();
	double getAverageLsr();
	void ReadData(); //czytanie danych
	void normal_vectors(); // wypelnie wektorow danymi, parzyste do lewego, nieparzyste do prawego
	void minus_vectors(); // wektory wypelniane obliczonymi danymi poprzez odjecie wartosci poprzedniego sampla od obecnego
	double first_calculation(double a, vector<INT16> b); // obliczanie tej sumy
	double first_calculation_minus(double a, vector<double> b); // funkcje x_minus to te same funkcje co wyzej tylko przymujace inny typ danych
	double entro(vector<INT16> a); // obliczanie entropi
	double entro_minus(vector<double> a); // Funkcja dokonuje rozk³adu LU macierzy A
	bool ludist(int n, double ** A); // Funkcja wyznacza wektor X na podstawie A i B
	bool lusolve(int n, double ** A, double * B, double * X);
	vector<double> counterRepeat(vector<INT16>canal, vector<double>vectorEPS);
	double SystemOfEquations(vector<INT16>canal);
	bool sign(double a);
	double EntroBit(vector<INT16>canal);
	double divideEPS(vector<INT16>canal);
};