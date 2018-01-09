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
	DWORD pcm; // WAVE type format
	WORD chanel; // mono/stereo
	DWORD sample_rate; // sample rate
	DWORD bytes_per_sec; // bytes/sec
	WORD block_alignment; // Block alignment
	WORD bits_per_sample; // Bits/sample
	FOURCC data; // "data" description header
	DWORD size_of_data; // size of data chunk
	const double eps = 1e-12; //wspolczynnik
	INT16 data_of_file;
	vector<INT16> leftCanal;
	vector<double> leftDifferential;
	vector<INT16> rightCanal;
	vector<double> rightDifferential;
	vector<double>minLsrVector;
	ofstream new_file;
	FILE *wf;
	double ammount_of_samples;
	double averageEPS;
	double averageBit;
	double averageLsr;
	int r;

	vector<double>predictCoderLeft;
	vector<double>predictCoderRight;

public:
	Wave(string name_of_wave, int r);
	~Wave();
	double getAverageEPS();
	double getAverageBit();
	double getAverageLsr();

	void ReadData(); 
	void NormalCanal(); 
	void DifferentialCoder(); 
	vector<double> predictCoder(vector<INT16>canal, vector<double>vectorEPS);

	double Energy(double a, vector<INT16> b); 
	double EnergyDifferential(double a, vector<double> b); 
	double Entropy(vector<INT16> a); 
	double EntropyDifferential(vector<double> a); 
	bool ludist(int n, double ** A); 
	bool lusolve(int n, double ** A, double * B, double * X);
	
	double SystemOfEquations(vector<INT16>canal);
	bool sign(double a);
	double EntroBit(vector<INT16>canal);
	double DivideEPS(vector<INT16>canal);
	void DecoderDifferential(vector<double>canal);
	void DecoderPredictive(vector<INT16>canal);
	void BothCanals();
};