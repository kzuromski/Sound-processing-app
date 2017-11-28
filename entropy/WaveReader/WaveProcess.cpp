#include "stdafx.h"
#include "WaveProcess.h"


WaveProcess::WaveProcess(string name_of_wave)
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
	ammount_of_samples = ((size_of_file - 48) - 2) / 2; // liczba probek z obu kanalow
	new_file << name_of_wave << endl;
	new_file << "Liczba próbek z obu kana³ów: " << fixed << ammount_of_samples << endl;
	FillingUpChanels(); // wektory wypelniane danymi
	double d = ((CalculateEnergy(left_chanel , ammount_of_samples / 2) + CalculateEnergy(right_chanel, ammount_of_samples / 2)) / 2);
	new_file << "Przeciêtna energia sygna³u: " << fixed << d << endl; // pierwsze obliczenia
	SamplesAfterScanning(); // wektory wypelniane obliczonymi danymi poprzez odjecie wartosci poprzedniej probki od obecnej
	d = ((CalculateEnergy(left_after_scan,ammount_of_samples / 2) + CalculateEnergy(right_after_scan,ammount_of_samples / 2)) / 2);
	new_file << "Przeciêtna energia sygna³u po skanowaniu ró¿nicowym: " << fixed << d << endl;  // drugie obliczenia
	new_file << "Entropia: " << ((CalculateEntrophy(left_chanel) + CalculateEntrophy(right_chanel)) / 2) << endl; // entro dla normalnych
	new_file << "Entropia z danych po skanowaniu ró¿nicowym: " << ((CalculateEntrophy(left_after_scan) + CalculateEntrophy(right_after_scan)) / 2) << endl; // entro dla tych po skalowaniu
	SystemOfEquations();
	
	average = CalculateEntrophy(eN);
	//new_file << "Entropia ze wspó³czynnikiem " << average << endl;
	
	//EntroBit();
}


WaveProcess::~WaveProcess()
{
	new_file.close();
	left_chanel.clear();
	right_chanel.clear();
	left_after_scan.clear();
	right_after_scan.clear();
	xDaszekVector.clear();
	eN.clear();
	maleA.clear();
}

void WaveProcess::ReadData()
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

void WaveProcess::FillingUpChanels()
{
	for (int i = 0; i < ammount_of_samples; i++)
	{
		fread(&data_of_file, sizeof(INT16), 1, wf);
		if (i % 2 == 0)
		{
			left_chanel.push_back(data_of_file);
		}
		else
		{
			right_chanel.push_back(data_of_file);
		}
	}
}

void WaveProcess::SamplesAfterScanning()
{
	for (int i = 0; i < left_chanel.size(); i++)
	{
		if (i == 0)
		{
			left_after_scan.push_back(left_chanel.at(i));
		}
		else
		{
			left_after_scan.push_back(left_chanel.at(i) - left_chanel.at(i - 1));
		}
	}
	for (int i = 0; i < right_chanel.size(); i++)
	{
		if (i == 0)
		{
			right_after_scan.push_back(right_chanel.at(i));
		}
		else
		{
			right_after_scan.push_back(right_chanel.at(i) - right_chanel.at(i - 1));
		}
	}
}

double WaveProcess::CalculateEnergy(vector<INT16> samples, double ammount_of_samples)
{
	double energy = 0;
	for (int i = 0; i < ammount_of_samples; i++)
	{
		energy += pow(samples.at(i),2);
	}

	energy = energy / ammount_of_samples;

	return energy;
}

double WaveProcess::CalculateEnergy(vector<double> samples, double ammount_of_samples)
{
	double energy = 0;
	for (int i = 0; i < ammount_of_samples; i++)
	{
		energy += pow(samples.at(i), 2);
	}

	energy = energy / ammount_of_samples;

	return energy;
}

double WaveProcess::CalculateEntrophy(vector<INT16> samples)
{
	double entrophy = 0;
	vector <double> buffor(131072, 0); // vector przetrzymujacy 2^17 miejsc
	for (int i = 0; i < samples.size(); i++)
	{
		buffor.at(samples.at(i) + 65536)++; // dodawanie powtarzajacych sie wartosci
	}
	for (int i = 0; i < 131072; i++)
	{
		if (buffor.at(i) != 0) // jezeli probka sie nie pojawila i jest 0 to nie mozemy jej uzyc, bo logarytm wywali nand
		{
			double p_i = (double)buffor.at(i) / samples.size(); // p_i = probka w indeksie 
			entrophy = entrophy + (p_i*log2(p_i));
		}
	}
	entrophy = entrophy * (-1);

	return entrophy;
}

double WaveProcess::CalculateEntrophy(vector<double> samples)
{
	double entrophy = 0;
	vector <double> buffor(131072, 0); // vector przetrzymujacy 2^17 miejsc
	for (int i = 0; i < samples.size(); i++)
	{
		buffor.at(samples.at(i) + 65536)++; // dodawanie powtarzajacych sie wartosci
	}
	for (int i = 0; i < 131072; i++)
	{
		if (buffor.at(i) != 0) // jezeli probka sie nie pojawila i jest 0 to nie mozemy jej uzyc, bo logarytm wywali nand
		{
			
			double p_i = (double)buffor.at(i) / samples.size(); // p_i = probka w indeksie i
			entrophy = entrophy + (p_i*log2(p_i));
		}
	}
	entrophy = entrophy * (-1);

	return entrophy;
}

bool WaveProcess::Ludist(int n, double ** A)
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

bool WaveProcess::Lusolve(int n, double ** A, double * B, double * X)
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

void WaveProcess::SystemOfEquations() {
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

	int N = ammount_of_samples / 2; //dla jednego kanalu
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
				sum += (right_chanel.at(z - i) * right_chanel.at(z - j) + left_chanel.at(z - i) * left_chanel.at(z - j)) * 0.5; //suma dla elementów macierzy X
				sum2 += (right_chanel.at(z) * right_chanel.at(z - i) + left_chanel.at(z) * left_chanel.at(z - i)) * 0.5; //suma dla elementów macierzy P
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

	int licznik = 0;
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			A[i][j] = a.at(licznik); // Macierz X (3x3)
			licznik++;
			B[i] = b.at(i); // Macierz P (3x1)
		}
	}

	//Warunek czy macierz ma wyznanik niezerowy
	if (Ludist(n, A) && Lusolve(n, A, B, X)) {}
	else cout << "DZIELNIK ZERO\n";


	for (int i = 0; i < r; i++) {
		maleA.push_back(X[i]);
	}

	xDaszekVector = NewXPushingAndCutting();

	for (int i = 0; i<ammount_of_samples / 2; i++)
	{
		eN.push_back((right_chanel.at(i) + left_chanel.at(i)) / 2 - xDaszekVector.at(i));
	}

	// usuwamy macierze z pamiêci
	for (i = 0; i < n; i++) delete[] A[i];
	delete[] A;
	delete[] B;
	delete[] X;
}

vector <double> WaveProcess::NewXPushingAndCutting()
{
	vector<double> xDaszekVector;
	double xDaszek = 0;
	for (size_t i = 0; i < ammount_of_samples / 2; i++)
	{
		xDaszek = 0;
		if (i == 0) {
			xDaszekVector.push_back((right_chanel.at(i) + left_chanel.at(i)) * 0.5);
		}
		if (i <= r)
		{
			xDaszekVector.push_back((right_chanel.at(i) + left_chanel.at(i)) * 0.5);
		}
		else {
			for (size_t j = 0; j < r; j++)
			{
				xDaszek += maleA.at(j) * ((right_chanel.at(i - j) + left_chanel.at(i - j)) * 0.5);
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
	return xDaszekVector;
}

double WaveProcess::GetAverage()
{
	return average;
}
