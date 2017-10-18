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
#include <map>
#include <algorithm>

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
	vector<INT16> left_minus;
	vector<INT16> right;
	vector<INT16> right_minus;
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
		ammount_of_samples = size_of_data / 4;

		new_file << name_of_wave << endl;
		new_file << "Size of data: " << size_of_data << endl;
		new_file << "Ammount of samples: " << ammount_of_samples << endl;
		normal_vectors(); // wektory wypelniane danymi
		new_file << "First calculation: " << ((first_calculation(ammount_of_samples / 2, left)+ first_calculation(ammount_of_samples / 2, right))/2) << endl; // pierwsze obliczenia

		minus_vectors(); // wektory wypelniane obliczonymi danymi poprzez odjecie wartosci poprzedniego sampla od obecnego
		new_file << "Second calculation: " << ((first_calculation(ammount_of_samples / 2, left_minus)+ first_calculation(ammount_of_samples / 2, right_minus))/2) << endl;  // drugie obliczenia

		
		new_file << "Third calculation: " << ((first_entro(left) + first_entro(right)) / 2) << endl; // entro dla normalnych

		new_file << "Fourth calculation: " << ((first_entro(left_minus) + first_entro(right_minus)) / 2) << endl; // entro dla tych errorow


	}
	 
	~WaveReader()
	{
		new_file.close();
	}

	void ReadData() //czytanie danych
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

	
	void normal_vectors() // wypelnie wektorow danymi
	{
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
	void minus_vectors() // wektory wypelniane obliczonymi danymi poprzez odjecie wartosci poprzedniego sampla od obecnego
	{
		for (int i = 0; i < left.size(); i++)
		{
			if (i == 0)
			{
				left_minus.push_back (left.at(i));
			}
			else
			{
				left_minus.push_back (left.at(i) - left.at(i - 1));
			}
		}
		for (int i = 0; i < right.size(); i++)
		{
			if (i == 0)
			{
				right_minus.push_back(right.at(i));
			}
			else
			{
				right_minus.push_back(right.at(i) - right.at(i - 1));
			}
		}
	}
	double first_calculation(double a, vector<INT16> b) // obliczanie tej sumy
	{
		double full = 0;
		for (double i = 0; i < a; i++)
		{
			full = (double)(full + (b.at(i) * b.at(i)));
		}

		full = full / a;

		return full;
	}

	double geting_max_value(vector<INT16> a)
	{
		double max_val_vec = *max_element(a.begin(), a.end());
		return max_val_vec;
	}
	double geting_min_value(vector<INT16> a)
	{
		double min_val_vec = *min_element(a.begin(), a.end());
		return min_val_vec;
	}

	double first_entro(vector<INT16> a)
	{
		double entro=0; // wartosc entro poczatkowo ustawiona na 0
		double min = geting_min_value(a); // najmniejsza wartosc w vektorze
		double max = geting_min_value(a); // najwieksza wartosc w wektorze

		for (min; min <= max; min++) // petla idaca od najmniejszego elementu do najwiekszego elementu
		{
			double same=0; // licznik dla takich powtorzen liczby
			for (int i = 0; i < a.size(); i++) // petla idaca od 0 do wielkosci vektora w ktorym sprawdzamy powtarzajace sie liczby
			{
				if (min == a.at(i)) // jezeli liczba ma wartosc min, czyli taka ktora teraz sprawdzamy to same ma sie powiekszyc
				{
					same++;
				}
			}
			double p_i = same / a.size(); // tutaj jest dzielenie ilosci powtorzen przez ilosc wszystkich elentow
			entro = entro + (p_i*log2(p_i)); // wzor entro czyli suma_elementow(pi*pi*log2pi)
		}

		return entro * (-1); // przed ta suma we wzorze byl jeszcze -
	}


	/*std::map<int, unsigned int> counter(const std::vector<INT16>& vals) {
		std::map<int, unsigned int> rv;

		for (auto val = vals.begin(); val != vals.end(); ++val) {
			rv[*val]++;
		}

		return rv;
	}

	void display(const std::map<int, unsigned int>& counts) {
		for (auto count = counts.begin(); count != counts.end(); ++count) {
			std::cout << "Value " << count->first << " has count "
				<< count->second << std::endl;
		}
	}*/



};

void main()
{
	WaveReader wave("sample.wav");
}