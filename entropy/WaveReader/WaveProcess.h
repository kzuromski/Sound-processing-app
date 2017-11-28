#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <Windows.h>
#include <fstream>
#include <math.h>
#include <iomanip>


using namespace std;

extern int r;

class WaveProcess
{
private:
	FOURCC riff; // "RIFF" file descripton header
	FOURCC wave; // "WAVE" file descripton header
	FOURCC fmt; // "fmt" description header
	FOURCC data; // "data" description header
	DWORD size_of_file; // size of file
	DWORD chunk; // size of WAVE section chunck
	DWORD sample_rate; // sample rate
	DWORD bytes_per_sec; // bytes/sec
	DWORD size_of_data; // size of data chunk
	WORD pcm; // WAVE type format
	WORD chanel; // mono/stereo
	WORD block_alignment; // Block alignment
	WORD bits_per_sample; // Bits/sample
	
	INT16 data_of_file;
	ofstream new_file;
	FILE *wf;
	double ammount_of_samples;
	const double eps = 1e-12; //wspolczynnik
	double average;

	vector<INT16> left_chanel;
	vector<INT16> right_chanel;
	vector<double> left_after_scan;
	vector<double> right_after_scan;
	vector<double> xDaszekVector;
	vector<double> eN;
	vector<double>maleA;
	
public:
	WaveProcess(string name_of_wave);
	~WaveProcess();
	void ReadData(); // czytanie danych z pliku
	void FillingUpChanels(); // wypelnianie kanalow danymi z pliku
	void SamplesAfterScanning();  // probki po skanowaniu roznicowym
	double CalculateEnergy(vector<INT16> samples, double ammount_of_samples); // obliczanie energii
	double CalculateEnergy(vector<double> samples, double ammount_of_samples); // obliczanie energii
	double CalculateEntrophy(vector<INT16> samples); // obliczanie entropi
	double CalculateEntrophy(vector<double> samples); // obliczanie entropi
	bool Ludist(int n, double ** A); // Funkcja dokonuje rozk³adu LU macierzy A
	bool Lusolve(int n, double ** A, double * B, double * X); // Funkcja wyznacza wektor X na podstawie A i B
	void SystemOfEquations(); // Funkcja obliczajaca 
	vector <double> NewXPushingAndCutting();
	double GetAverage();
	
	/*void EntroBit() {

	for (int i = 0; i < maleA.size(); i++) {
	cout << maleA.at(i) << endl;
	}
	auto max = max_element(begin(maleA), end(maleA));
	auto min = min_element(begin(maleA), end(maleA));

	if (abs(*min) > *max)
	*max = *min;

	auto position = distance(begin(maleA), max);
	cout << "max " << *max << endl;
	cout << "position " << position << endl;
	cout << endl;
	}*/

};

