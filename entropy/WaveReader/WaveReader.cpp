#include "stdafx.h"
#include "windows.h"
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
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
	vector<double> left_minus;
	vector<INT16> right;
	vector<double> right_minus;
	ofstream new_file;
	FILE *wf;
	double ammount_of_samples;
	const double eps = 1e-12; //wspolczynnik

	vector <double> xDaszekVector;
	vector <double> eN;
	vector<double>maleA;
	double average;
	int r = 3;

public:
	WaveReader(string name_of_wave)
	{
		cout << name_of_wave << endl;
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
		new_file << "Liczba próbek z obu kana³ów: " << fixed << ammount_of_samples << endl;
		normal_vectors(); // wektory wypelniane danymi

		//double d = ((first_calculation(ammount_of_samples / 2, left) + first_calculation(ammount_of_samples / 2, right)) / 2);
		//new_file << "Przeciêtna energia sygna³u: " << fixed << d << endl; // pierwsze obliczenia
		//minus_vectors(); // wektory wypelniane obliczonymi danymi poprzez odjecie wartosci poprzedniej probki od obecnej
		//d = ((first_calculation_minus(ammount_of_samples / 2, left_minus) + first_calculation_minus(ammount_of_samples / 2, right_minus)) / 2);
		//new_file << "Przeciêtna energia sygna³u po skanowaniu ró¿nicowym: " << fixed << d << endl;  // drugie obliczenia
		//
		//new_file << "Entropia: " << ((entro(left) + entro(right)) / 2) << endl; // entro dla normalnych
		//new_file << "Entropia z danych po skanowaniu ró¿nicowym: " << ((entro_minus(left_minus) + entro_minus(right_minus))/2) << endl; // entro dla tych po skalowaniu
		SystemOfEquations();
		average = entro_minus(eN);
		new_file << "Entropia ze wspó³czynnikiem " << average <<endl;
		EntroBit();

	}
	 
	~WaveReader()
	{
		new_file.close();
	}

	double getAverage() {
		return average;
	}

	int getR() {
		return r;
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
		for (INT32 i = 0; i < a; i++)
		{
			full = (double)(full + (((double)b.at(i) * (double)b.at(i))));
		}

		full = full / a;
		
		return full;
	}

	double first_calculation_minus(double a, vector<double> b) // funkcje x_minus to te same funkcje co wyzej tylko przymujace inny typ danych
	{
		double full = 0;
		for (INT32 i = 0; i < a; i++)
		{
			full = (double)(full + (b.at(i) * b.at(i)));
		}

		full = full / a;

		return full;
	}

	double entro(vector<INT16> a) // obliczanie entropi
	{
		double entro = 0;
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
				entro = entro + (p_i*log2(p_i));
			}
		}
		return entro *(-1);
	}

	double entro_minus(vector<double> a)
	{
		double entro = 0;
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
				entro = entro + (p_i*log2(p_i));
			}
		}
		return entro *(-1);
	}

	// Funkcja dokonuje rozk³adu LU macierzy A
	bool ludist(int n, double ** A)
	{
		int i, j, k;

		for (k = 0; k < n - 1; k++)
		{
			if (fabs(A[k][k]) < eps) return false;

			for (i = k + 1; i < n; i++)
				A[i][k] /= A[k][k];

			for (i = k + 1; i < n; i++)
				for (j = k + 1; j < n; j++)
					A[i][j] -= A[i][k] * A[k][j];
		}

		return true;
	}

	// Funkcja wyznacza wektor X na podstawie A i B
	bool lusolve(int n, double ** A, double * B, double * X)
	{
		int    i, j;
		double s;

		X[0] = B[0];

		for (i = 1; i < n; i++)
		{
			s = 0;

			for (j = 0; j < i; j++) s += A[i][j] * X[j];

			X[i] = B[i] - s;
		}

		if (fabs(A[n - 1][n - 1]) < eps) return false;

		X[n - 1] /= A[n - 1][n - 1];

		for (i = n - 2; i >= 0; i--)
		{
			s = 0;

			for (j = i + 1; j < n; j++) s += A[i][j] * X[j];

			if (fabs(A[i][i]) < eps) return false;

			X[i] = (X[i] - s) / A[i][i];
		}

		return true;
	}

	
	void SystemOfEquations() {
		double **A, *B, *X;
		int n, i, j;
		
		n = r;

		cout << setprecision(25) << fixed;

		// tworzymy macierze A, B i X
		A = new double *[n];
		B = new double[n];
		X = new double[n];

		for (i = 0; i < n; i++) 
			A[i] = new double[n];

		int N = ammount_of_samples/2; //dla jednego kanalu
		double sum = 0; //zerowanie sumy
		double sum2 = 0;
		vector<double>a; //wektor dla macierzy X
		vector<double>b; //wektor dla macierzy P

		for (i = 1; i <= r; i++)
		{
			for (j = 1; j <= r; j++)
			{
				for (int z = r; z < N; z++)
				{
					sum += (right.at(z - i) * right.at(z - j) + left.at(z - i) * left.at(z - j)) * 0.5; //suma dla elementów macierzy X
					sum2 += (right.at(z) * right.at(z - i) + left.at(z) * left.at(z - i)) * 0.5; //suma dla elementów macierzy P
				}
				a.push_back(sum);

				if (j == 1)
				{
					b.push_back(sum2);
				}
				sum = 0;
				sum2 = 0;
			}
		}

		int licznik=0;
		for (i = 0; i < n; i++){
			for (j = 0; j < n; j++) {
				A[i][j] = a.at(licznik); // Macierz X (3x3)
				licznik++;
				B[i] = b.at(i); // Macierz P (3x1)
			}
		}

		//Warunek czy macierz ma wyznanik niezerowy
		if (ludist(n, A) && lusolve(n, A, B, X)){}
		else cout << "DZIELNIK ZERO\n";

		
		for (int i = 0; i < r; i++) {
			maleA.push_back(X[i]);
		}

		double xDaszek=0;
		for (size_t i = 0; i < ammount_of_samples/2; i++)
		{
			xDaszek = 0;
			if (i == 0) {
				xDaszekVector.push_back((right.at(i) + left.at(i)) * 0.5);
			}
			else if (i <= r)
			{
				xDaszekVector.push_back((right.at(i - 1) + left.at(i - 1)) * 0.5);
			}
			else {
				for (size_t j = 0; j < r; j++)
				{
					xDaszek += maleA.at(j) * ((right.at(i - j) + left.at(i - j)) * 0.5);
				}
				if (xDaszek > pow(2, 15) - 1)
				{
					xDaszek = pow(2, 15) - 1;
				}
				else if (xDaszek < -pow(2, 15))
				{
					xDaszek = -pow(2, 15);
				}
				xDaszekVector.push_back(floor(xDaszek + 0.5));
			}
		}

		for (int i = 0; i<ammount_of_samples / 2; i++)
		{
			eN.push_back((right.at(i)+left.at(i))/2 - xDaszekVector.at(i));
		}

		// usuwamy macierze z pamiêci
		for (i = 0; i < n; i++) delete[] A[i];
		delete[] A;
		delete[] B;
		delete[] X;
	}

	bool sign(double a){
		if (a >= 0)
			return 1;
		else
			return 0;
	}

	void EntroBit() {

		auto max = max_element(begin(maleA), end(maleA));
		auto min = min_element(begin(maleA), end(maleA));
		auto position = 0; 

		bool positive = true;
		if (abs(*min) > *max) {
			position = distance(begin(maleA), min);
			*max = abs(*min); 
			positive = false;
		}
		else {
			position = distance(begin(maleA), max);
		}

		*max = float(*max); 
		vector<int>si;
		vector<double>decod;
		vector<double>aDaszek;
		vector<double>xDaszekVector2;
		vector<double>eN2;
		
		double Lsr;
		double minLsr = 100;
		int diagramBit = 0;
		for (int b = 8; b <= 24; b++) {

			for (int i = 0; i < maleA.size(); i++) {
				aDaszek.push_back(floor(abs(maleA.at(i)) / (*max) * (pow(2, b) - 1) + 1 / 2));//liczy dobrze
				si.push_back(sign(maleA.at(i)));
			}

			for (int i = 0; i < maleA.size(); i++) {
				decod.push_back((aDaszek.at(i) / (pow(2, b) - 1) * (*max)) * (si.at(i) * 2 - 1));//liczby dobrze
			}

			double xDaszek2 = 0;
			for (size_t i = 0; i < ammount_of_samples / 2; i++)
			{
				xDaszek2 = 0;
				if (i == 0) {
					xDaszekVector2.push_back((right.at(i) + left.at(i)) * 0.5);
				}
				else if (i <= r)
				{
					xDaszekVector2.push_back((right.at(i - 1) + left.at(i - 1)) * 0.5);
				}
				else {
					for (size_t j = 0; j < r; j++)
					{
						xDaszek2 += decod.at(j) * ((right.at(i - j) + left.at(i - j)) * 0.5);
					}
					if (xDaszek2 > pow(2, 15) - 1)
					{
						xDaszek2 = pow(2, 15) - 1;
					}
					else if (xDaszek2 < -pow(2, 15))
					{
						xDaszek2 = -pow(2, 15);
					}
					xDaszekVector2.push_back(floor(xDaszek2 + 0.5));
				}
			}

			for (int i = 0; i<ammount_of_samples / 2; i++)
			{
				eN2.push_back((right.at(i) + left.at(i)) / 2 - xDaszekVector2.at(i));
			}

			Lsr = entro_minus(eN2) + ((32 + (r - 1) * (b + 1) + 10) / ammount_of_samples);//czasami rosnie, czasami nie xD
			if (minLsr > Lsr) {
				minLsr = Lsr;
				diagramBit = b;
			}

			si.clear();
			eN2.clear();
			decod.clear();
			aDaszek.clear();
			xDaszekVector2.clear();
		}	
		cout << "bit: " << diagramBit << endl;
	}
};

void main()
{
	WaveReader wave1("ATrain.wav");
	WaveReader wave2("BeautySlept.wav");
	WaveReader wave3("death2.wav");
	WaveReader wave4("experiencia.wav");
	WaveReader wave5("chanchan.wav");
	WaveReader wave6("female_speech.wav");
	WaveReader wave7("FloorEssence.wav");
	WaveReader wave8("ItCouldBeSweet.wav");
	WaveReader wave9("Layla.wav");
	WaveReader wave10("LifeShatters.wav");
	WaveReader wave11("macabre.wav");
	WaveReader wave12("male_speech.wav");
	WaveReader wave13("SinceAlways.wav");
	WaveReader wave14("thear1.wav");
	WaveReader wave15("TomsDiner.wav");
	WaveReader wave16("velvet.wav");

	/*double tab[16];
	double average = 0;
	
	tab[0] = wave1.getAverage();
	tab[1] = wave2.getAverage();
	tab[2] = wave3.getAverage();
	tab[3] = wave4.getAverage();
	tab[4] = wave5.getAverage();
	tab[5] = wave6.getAverage();
	tab[6] = wave7.getAverage();
	tab[7] = wave8.getAverage();
	tab[8] = wave9.getAverage();
	tab[9] = wave10.getAverage();
	tab[10] = wave11.getAverage();
	tab[11] = wave12.getAverage();
	tab[12] = wave13.getAverage();
	tab[13] = wave14.getAverage();
	tab[14] = wave15.getAverage();
	tab[15] = wave16.getAverage();
	
	for (int i = 0; i < 16; i++) {
		average += tab[i];
	}
	average /= 16;*/

	//ofstream averageFile("average.txt", ios::app);
	//averageFile << wave1.getR() << " " << average << endl;
	system("pause");
}