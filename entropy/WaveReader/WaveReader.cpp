// WaveReader.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <cstdio>
#include <cmath>
#include <string>
#include "windows.h"
#include "mmsystem.h"
#include <fstream>
#include <time.h>
#include <vector>

using namespace std;

class WaveReader
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
	INT16 data_of_file;
	vector<INT16> left;
	vector<INT16> right;
	ofstream new_file;
	FILE *wf;
	double ammount_of_samples;
public:
	WaveReader(string name_of_wave)
	{
		
		const char * c = name_of_wave.c_str();
		wf = fopen(c, "r");

		if (wf == NULL)
		{
			cout << "File is not opened" << endl;
			exit(-1);
		}

		SYSTEMTIME time;
		GetLocalTime(&time);
		string file_name;
		file_name = to_string(time.wHour) + "_" + to_string(time.wMinute) + "_" + to_string(time.wSecond);
		ReadData();

		ofstream new_file(file_name + ".txt");
		double ammount_of_samples = size_of_data / 4;
		new_file << name_of_wave << endl;
		new_file << "Size of data: " << size_of_data << endl;
		new_file << "Ammount of samples: " << ammount_of_samples << endl;
		calculations();
		new_file << "First calculation of left side: " << first_calculation(ammount_of_samples / 2, left) << endl;
		new_file << "First calculation of right side: " << first_calculation(ammount_of_samples / 2, right) << endl;

	}
	 
	~WaveReader()
	{
		new_file.close();
	}

	void ReadData() //czytanie wartoœci próbek do wektorów kana³ów
	{
		fread(&riff, sizeof(FOURCC), 1, wf);
		fread(&size_of_file, sizeof(DWORD), 1, wf);
		fread(&wave, sizeof(FOURCC), 1, wf);
		fread(&fmt, sizeof(FOURCC), 1, wf);
		fread(&chunk, sizeof(FOURCC), 1, wf);
		fread(&pcm, sizeof(WORD), 1, wf);
		fread(&chanel, sizeof(WORD), 1, wf);
		fread(&sample_rate, sizeof(DWORD), 1, wf);
		fread(&bytes_per_sec, sizeof(DWORD), 1, wf);
		fread(&block_alignment, sizeof(WORD), 1, wf);
		fread(&bits_per_sample, sizeof(WORD), 1, wf);
		fread(&data, sizeof(FOURCC), 1, wf);
		fread(&size_of_data, sizeof(DWORD), 1, wf);
	}

	
	void calculations()
	{
		ammount_of_samples = size_of_data / 4;
		for (int i = 0; i < ammount_of_samples; i++)
		{
			fread(&data_of_file, sizeof(INT16), 1, wf);
			if (i % 2 == 0)
			{
				left.push_back(data_of_file);
			}
			else
			{
				right.push_back(data_of_file);
			}
		}

		
	}

	double first_calculation(double a, vector<INT16> b)
	{
		double full = 0;
		for (double i = 0; i < a; i++)
		{
			full = (double)(full + (b.at(i) * b.at(i)));
		}

		full = full / a;

		return full;
	}

};

void main()
{
	WaveReader wave("gun.wav");
}