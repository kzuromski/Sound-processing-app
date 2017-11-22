// WaveReader.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <cstdio>
#include <cmath>
#include <string>
#include "windows.h"
#include <fstream>
#include <vector>
#include <algorithm>
#include <iomanip>

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
	vector< vector<INT16> >canals;
	vector< vector<double> >canalsMinus;
	vector< vector<double> > xDaszekVector;
	double xDaszek;
	vector <double> errorN;
	ofstream new_file;
	FILE *wf;
	double ammount_of_samples;
	const double eps = 1e-12; //wspolczynnik

public:
	WaveReader(string name_of_wave)
	{
		clock_t start = clock();
		canals.push_back(vector<INT16>());
		canals.push_back(vector<INT16>());
		canalsMinus.push_back(vector<double>());
		canalsMinus.push_back(vector<double>());
		const char * c = name_of_wave.c_str();
		wf = fopen(c, "rb");

		if (wf == NULL)
		{
			cout << "File is not opened" << endl;
			exit(-1);
		}

		ReadData();

		ofstream new_file(name_of_wave + ".txt");
		ammount_of_samples = ((size_of_file - 48)-2)/ 2; // liczba probek z obu kanalow
		new_file << name_of_wave << endl;
		cout << name_of_wave << endl;
		new_file << "Liczba próbek z obu kana³ów: " << fixed << ammount_of_samples << endl;
		normal_vectors(); // wektory wypelniane danymi

		double d = ((first_calculation(ammount_of_samples / 2, canals[0]) + first_calculation(ammount_of_samples / 2, canals[1])) / 2);
		new_file << "Przeciêtna energia sygna³u: " << fixed << d << endl; // pierwsze obliczenia
		minus_vectors(); // wektory wypelniane obliczonymi danymi poprzez odjecie wartosci poprzedniej probki od obecnej
		d = ((first_calculation(ammount_of_samples / 2, canalsMinus[0]) + first_calculation(ammount_of_samples / 2, canalsMinus[1])) / 2);
		new_file << "Przeciêtna energia sygna³u po skanowaniu ró¿nicowym: " << fixed << d << endl;  // drugie obliczenia
		
		new_file << "Entropia: " << ((entropy(canals[0]) + entropy(canals[1])) / 2) << endl; // entro dla normalnych
		new_file << "Entropia z danych po skanowaniu ró¿nicowym: " << ((entropy(canalsMinus[0]) + entropy(canalsMinus[1])) / 2) << endl; // entro dla tych po skalowaniu
		SystemOfEquations();
		new_file << "Entropia ze wspó³czynnikiem " << entropy(errorN)<<endl;
		cout << "Czas dziaania programu dla pliku(w sekundach): " << (clock() - start) / CLOCKS_PER_SEC << endl;
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

	
	void normal_vectors() // wypelnie wektorow danymi, parzyste do lewego, nieparzyste do prawego
	{
		for (int i = 0; i < ammount_of_samples; i++)
		{
			fread(&data_of_file, sizeof(INT16), 1, wf);
			if (i % 2 == 0)
			{
				canals[0].push_back(data_of_file);
			}
			else
			{
				canals[1].push_back(data_of_file);
			}
		}
	}

	void minus_vectors() // wektory wypelniane obliczonymi danymi poprzez odjecie wartosci poprzedniego sampla od obecnego
	{
		for (int i = 0; i < canals[0].size(); i++)
		{
			if (i == 0)
			{
				canalsMinus[0].push_back (canals[0].at(i));
			}
			else
			{
				canalsMinus[0].push_back (canals[0].at(i) - canals[0].at(i - 1));
			}
		}
		for (int i = 0; i < canals[1].size(); i++)
		{
			if (i == 0)
			{
				canalsMinus[1].push_back(canals[1].at(i));
			}
			else
			{
				canalsMinus[1].push_back(canals[1].at(i) - canals[1].at(i - 1));
			}
		}
	}

	double first_calculation(double a, vector<INT16> b) // obliczanie tej sumy
	{
		double full = 0;
		for (INT32 i = 0; i < a; i++)
		{
			full = (double)(full + (((double)b.at(i) * (double)b.at(i))));
		}

		full = full / a;
		
		return full;
	}

	double first_calculation(double a, vector<double> b) // funkcje x_minus to te same funkcje co wyzej tylko przymujace inny typ danych
	{
		double full = 0;
		for (INT32 i = 0; i < a; i++)
		{
			full = (double)(full + (b.at(i) * b.at(i)));
		}

		full = full / a;

		return full;
	}

	double entropy(vector<INT16> a) // obliczanie entropi
	{
		double entropy = 0;
		vector <double> buffor(131072, 0); // vector przetrzymujacy 2^17 miejsc
		for (int i = 0; i < a.size(); i++)
		{
			buffor.at(a.at(i) + 65536)++; // dodawanie powtarzajacych sie wartosci
		}
		for (int i = 0; i < 131072; i++)
		{
			if (buffor.at(i) != 0) // jezeli probka sie nie pojawila i jest 0 to nie mozemy jej uzyc, bo logarytm wywali nand
			{
				double p_i = (double)buffor.at(i) / a.size();
				entropy = entropy + (p_i*log2(p_i));
			}
		}
		return entropy *(-1);
	}

	double entropy(vector<double> a)
	{
		double entropy = 0;
		vector <double> buffor(131072, 0);
		for (int i = 0; i < a.size(); i++)
		{
			buffor.at(a.at(i)+ 65536)++;
		}
		for (int i = 0; i < 131072; i++)
		{
			if (buffor.at(i) != 0)
			{
				double p_i = (double)buffor.at(i) / a.size();
				entropy = entropy + (p_i*log2(p_i));
			}
		}
		return entropy *(-1);
	}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



	// Funkcja dokonuje rozk³adu LU macierzy A
	bool ludist(int n, double ** A)
	{
		for (int k = 0; k < n - 1; k++)
		{
			if (fabs(A[k][k]) < eps)
			{
				return false;
			}

			for (int i = k + 1; i < n; i++)
			{
				A[i][k] /= A[k][k];
			}

			for (int i = k + 1; i < n; i++)
			{
				for (int j = k + 1; j < n; j++)
				{
					A[i][j] -= A[i][k] * A[k][j];
				}
			}
		}
		return true;
	}

	// Funkcja wyznacza wektor X na podstawie A i B
	bool lusolve(int n, double ** A, double * B, double * X)
	{
		int i, j;
		double s;

		X[0] = B[0];

		for (i = 1; i < n; i++)
		{
			s = 0;

			for (j = 0; j < i; j++)
			{
				s += A[i][j] * X[j];
			}
			X[i] = B[i] - s;
		}

		if (fabs(A[n - 1][n - 1]) < eps)
		{
			return false;
		}

		X[n - 1] /= A[n - 1][n - 1];

		for (i = n - 2; i >= 0; i--)
		{
			s = 0;

			for (j = i + 1; j < n; j++)
			{
				s += A[i][j] * X[j];
			}

			if (fabs(A[i][i]) < eps)
			{
				return false;
			}
			X[i] = (X[i] - s) / A[i][i];
		}
		return true;
	}

	vector<double> calculateError(vector<double> valuesofA, double * X, int canal, int r)
	{
		for (int i = 0; i < r; i++)
		{
			valuesofA.push_back(X[i]);
		}

		for (size_t i = 0; i < ammount_of_samples / 2; i++)
		{
			xDaszek = 0;
			if (i == 0)
			{
				xDaszekVector[canal].push_back(canals[canal].at(i));
			}
			else if (i <= r)
			{
				xDaszekVector[canal].push_back(canals[canal].at(i - 1));
			}
			else
			{
				for (size_t j = 0; j < r; j++)
				{
					xDaszek += valuesofA.at(j) * canals[canal].at(i - j);
				}
				if (xDaszek > pow(2, 15) - 1)
				{
					xDaszek = pow(2, 15) - 1;
				}
				else if (xDaszek < -pow(2, 15))
				{
					xDaszek = -pow(2, 15);
				}
				xDaszekVector[canal].push_back(floor(xDaszek + 0.5));
			}
		}

		for (int i = 0; i<ammount_of_samples / 2; i++)
		{
			errorN.push_back(canals[canal].at(i) - xDaszekVector[canal].at(i));
		}
		return errorN;
	}

	
	void SystemOfEquations()
	{
		int i, j;
		int r = 6;
		int n = r;
		double **A, *B, *X; // tworzymy macierze A, B i X
		A = new double *[n];
		B = new double[n];
		X = new double[n];
		double sumOfX = 0; //zerowanie sumy
		double sumOfP = 0;
		vector<double>valuesOfX; //wektor dla macierzy X
		vector<double>valuesOfP; //wektor dla macierzy P
		int N = ammount_of_samples / 2; //Iloœæ sampli dla jednego kana³u
		vector<double>valuesofA;
		int licznik = 0;
		
		cout << setprecision(4) << fixed;
		for (i = 0; i < n; i++)
		{
			A[i] = new double[n];
		}

		for (j = 1; j <= r; j++)
		{
			for (i = 1; i <= r; i++)
			{
				for (int z = r; z < N; z++)
				{
					sumOfX += canals[1].at(z - i) * canals[1].at(z - j); //Suma dla elementów macierzy X
					sumOfP += canals[1].at(z) * canals[1].at(z - i); //Suma dla elementów macierzy P
				}
				valuesOfX.push_back(sumOfX);
				valuesOfP.push_back(sumOfP);
				sumOfX = 0;
				sumOfP = 0;
			}
		}
		
		for (j = 0; j < n; j++)
		{
			for (i = 0; i < n; i++)
			{
				A[i][j] = valuesOfX.at(licznik); // Macierz X (3x3)
				licznik++;
				B[i] = valuesOfP.at(i); // Macierz P (3x1)
			}
		}
		

		// rozwi¹zujemy uk³ad i wyœwietlamy wyniki
		double sumA = 0;
		if (ludist(n, A) && lusolve(n, A, B, X))
		{
			for (i = 0; i < n; i++) 
			{
				sumA += X[i];
			}
		}
		else cout << "DZIELNIK ZERO\n";

		xDaszekVector.push_back(vector<double>());
		xDaszekVector.push_back(vector<double>());

		vector<double> errorN1 = calculateError(valuesofA, X, 0, r);
		vector<double> errorN2 = calculateError(valuesofA, X, 1, r);

		for (int i = 0; i < ammount_of_samples/2; i++)
		{
			errorN.at(i) = (errorN1.at(i) + errorN2.at(i))/2;
		}

		// usuwamy macierze z pamiêci
		for (i = 0; i < n; i++)
		{
			delete[] A[i];
		}
		delete[] A;
		delete[] B;
		delete[] X;
	}
};

void main()
{
	clock_t start = clock();
	cout << "Przetwarzanie plikow: " << endl;
	WaveReader wave1("ATrain.wav");
	WaveReader wave2("BeautySlept.wav");
	WaveReader wave3("death2.wav");
	WaveReader wave4("experiencia.wav");
	WaveReader wave5("female_speech.wav");
	WaveReader wave6("FloorEssence.wav");
	WaveReader wave7("ItCouldBeSweet.wav");
	WaveReader wave8("Layla.wav");
	WaveReader wave9("LifeShatters.wav");
	WaveReader wave10("macabre.wav");
	WaveReader wave11("male_speech.wav");
	WaveReader wave12("SinceAlways.wav");
	WaveReader wave13("thear1.wav");
	WaveReader wave14("TomsDiner.wav");
	WaveReader wave15("velvet.wav");
	WaveReader wave16("chanchan.wav");
	cout << "Czas dzia³ania programu(w sekundach): " << (clock() - start)/CLOCKS_PER_SEC << endl;
	system("pause");
}