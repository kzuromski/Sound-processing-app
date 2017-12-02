#include "Wave.h"

Wave::Wave(string name_of_wave, int r)
{
	this->r = r;
	cout << name_of_wave << endl;
	name_of_wave = ".\\Waves\\" + name_of_wave;
	const char *c = name_of_wave.c_str();

	wf = fopen(c, "rb");
	if (wf == NULL) {
		cout << "File is not opened" << endl;
		exit(-1);
	}
	ReadData();

	ofstream new_file(name_of_wave + ".txt");
	ammount_of_samples = ((size_of_file - 48) - 2) / 2; // liczba probek z obu kanalow
	new_file << name_of_wave << endl;
	new_file << "Liczba próbek z obu kana³ów: " << fixed << ammount_of_samples << endl;
	normal_vectors(); // wektory wypelniane danymi

	//double d = ((first_calculation(ammount_of_samples / 2, left) + first_calculation(ammount_of_samples / 2, right)) / 2);
	//new_file << "Przeciêtna energia sygna³u: " << fixed << d << endl; // pierwsze obliczenia
	//minus_vectors(); // wektory wypelniane obliczonymi danymi poprzez odjecie wartosci poprzedniej probki od obecnej
	//d = ((first_calculation_minus(ammount_of_samples / 2, left_minus) + first_calculation_minus(ammount_of_samples / 2, right_minus)) / 2);
	//new_file << "Przeciêtna energia sygna³u po skanowaniu ró¿nicowym: " << fixed << d << endl;  // drugie obliczenia
	//new_file << "Entropia: " << ((entro(left) + entro(right)) / 2) << endl; // entro dla normalnych
	//new_file << "Entropia z danych po skanowaniu ró¿nicowym: " << ((entro_minus(left_minus) + entro_minus(right_minus)) / 2) << endl; // entro dla tych po skalowaniu

	//r<2;40>
	//averageEPS = (SystemOfEquations(right) + SystemOfEquations(left)) / 2;
	//new_file << "Entropia ze wspó³czynnikiem " << averageEPS << endl;

	//r<10;600>
	//averageBit = ceil((EntroBit(right) + EntroBit(left)) / 2);

	//r<120;4>
	averageLsr = (divideEPS(right) + divideEPS(left)) / 2;
}

Wave::~Wave() {
	new_file.close();
}

double Wave::getAverageEPS() {
	return averageEPS;
}

double Wave::getAverageBit() {
	return averageBit;
}

double Wave::getAverageLsr() {
	return averageLsr;
}

void Wave::ReadData() 
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

void Wave::normal_vectors() 
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

void Wave::minus_vectors()
{
	for (int i = 0; i < left.size(); i++)
	{
		if (i == 0)
		{
			left_minus.push_back(left.at(i));
		}
		else
		{
			left_minus.push_back(left.at(i) - left.at(i - 1));
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

double Wave::first_calculation(double a, vector<INT16> b) 
{
	double full = 0;
	for (INT32 i = 0; i < a; i++)
	{
		full = (double)(full + (((double)b.at(i) * (double)b.at(i))));
	}

	full = full / a;

	return full;
}

double Wave::first_calculation_minus(double a, vector<double> b) 
{
	double full = 0;
	for (INT32 i = 0; i < a; i++)
	{
		full = (double)(full + (b.at(i) * b.at(i)));
	}

	full = full / a;

	return full;
}

double Wave::entro(vector<INT16> a) 
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

double Wave::entro_minus(vector<double> a)
{
	double entro = 0;
	vector <double> buffor(131072, 0);
	for (int i = 0; i < a.size(); i++)
	{
		buffor.at(a.at(i) + 65536)++;
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

bool Wave::ludist(int n, double ** A)
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

bool Wave::lusolve(int n, double ** A, double * B, double * X)
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

vector<double> Wave::sendEntropia(vector<INT16>canal, vector<double>vectorEPS) {
	double sumError = 0;
	vector <double> entropiaEPS;
	vector <double> absoluteError;
	for (size_t i = 0; i < ammount_of_samples / 2; i++) {
		sumError = 0;
		if (i == 0)
			absoluteError.push_back(canal.at(i));
		if (i <= r)
			absoluteError.push_back(canal.at(i));
		else {
			for (size_t j = 0; j < r; j++)
				sumError += vectorEPS.at(j) * canal.at(i - j);
			if (sumError > pow(2, 15) - 1)
				sumError = pow(2, 15) - 1;
			else if (sumError < -pow(2, 15))
				sumError = -pow(2, 15);
			absoluteError.push_back(floor(sumError + 0.5));
		}
	}
	for (int i = 0; i < ammount_of_samples / 2; i++)
		entropiaEPS.push_back(canal.at(i) - absoluteError.at(i));
	
	return entropiaEPS;
}

double Wave::SystemOfEquations(vector<INT16>canal) {
	double **A, *B, *X;
	int n, i, j;
	n = r;

	cout << setprecision(10) << fixed;
	A = new double *[n];
	B = new double[n];
	X = new double[n];

	for (i = 0; i < n; i++)
		A[i] = new double[n];

	int N = ammount_of_samples / 2;
	double sumX = 0;
	double sumP = 0;
	vector<double>matrixX;
	vector<double>matrixP;
	vector<double> vectorEPS;

	for (i = 1; i <= r; i++) {
		for (j = 1; j <= r; j++) {
			for (int z = r; z < N; z++) {
				sumX += canal.at(z - i) * canal.at(z - j);
				sumP += canal.at(z) * canal.at(z - i);
			}

			if (j == 1)
				matrixP.push_back(sumP);
			matrixX.push_back(sumX);
			sumX = 0;
			sumP = 0;
		}
	}

	int licznik = 0;
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			A[i][j] = matrixX.at(licznik);
			licznik++;
			B[i] = matrixP.at(i);
		}
	}

	if (ludist(n, A) && lusolve(n, A, B, X)) {}
	else cout << "DZIELNIK ZERO\n";

	for (int i = 0; i < r; i++) 
		vectorEPS.push_back(X[i]);

	vector<double> Entropia;
	Entropia = sendEntropia(canal, vectorEPS);
	double returnEntropia = entro_minus(Entropia);

	for (i = 0; i < n; i++)
		delete[] A[i];
	delete[] A;
	delete[] B;
	delete[] X;

	return returnEntropia;
}

bool Wave::sign(double a) {
	if (a >= 0)
		return 1;
	else
		return 0;
}

int Wave::EntroBit(vector<INT16>canal) {

	double **A, *B, *X;
	int n, i, j;
	n = r;
	int N = ammount_of_samples / 2;

	cout << setprecision(10) << fixed;
	A = new double *[n];
	B = new double[n];
	X = new double[n];

	for (i = 0; i < n; i++)
		A[i] = new double[n];

	double sumX = 0;
	double sumP = 0;
	vector<double>matrixX;
	vector<double>matrixP;
	vector<double>vectorEPS;

	for (i = 1; i <= r; i++) {
		for (j = 1; j <= r; j++) {
			for (int z = r ; z < N ; z++) {
				sumX += canal.at(z - i) * canal.at(z - j);
				sumP += canal.at(z) * canal.at(z - i);
			}

			if (j == 1)
				matrixP.push_back(sumP);
			matrixX.push_back(sumX);
			sumX = 0;
			sumP = 0;
		}
	}

	int licznik = 0;
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			A[i][j] = matrixX.at(licznik);
			licznik++;
			B[i] = matrixP.at(i);
		}
	}

	if (ludist(n, A) && lusolve(n, A, B, X)) {}
	else cout << "DZIELNIK ZERO\n";

	for (int i = 0; i < r; i++) 
		vectorEPS.push_back(X[i]);

	for (i = 0; i < n; i++)
		delete[] A[i];
	delete[] A;
	delete[] B;
	delete[] X;

	auto max = max_element(begin(vectorEPS), end(vectorEPS));
	auto min = min_element(begin(vectorEPS), end(vectorEPS));
	auto position = 0;

	bool positive = true;
	if (abs(*min) > *max) {
		position = distance(begin(vectorEPS), min);
		*max = abs(*min);
		positive = false;
	}
	else {
		position = distance(begin(vectorEPS), max);
	}

	*max = float(*max);
	vector<int>si;
	vector<double>coder;
	vector<double>decoder;

	double Lsr;
	double entropia = 0;
	double minLsr = 100;
	int diagramBit = 0;
	for (int b = 8; b <= 24; b++) {

		for (int i = 0; i < r; i++) {
			coder.push_back(floor(abs(vectorEPS.at(i)) / (*max) * (pow(2, b) - 1) + 0.5));
			si.push_back(sign(vectorEPS.at(i)));
		}

		for (int i = 0; i < r; i++)
			decoder.push_back(((coder.at(i) / (pow(2, b) - 1)) * (*max)) * (si.at(i) * 2 - 1));

		vector<double> decoderEntropia;
		decoderEntropia = sendEntropia(canal, decoder);

		Lsr = entro_minus(decoderEntropia) + ((32 + (r - 1) * (b + 1) + 10) / ammount_of_samples);
		if (minLsr > Lsr) {
			minLsr = Lsr;
			diagramBit = b;
		}	

		si.clear();
		decoder.clear();
		coder.clear();
	}
	
	return diagramBit;
}

double Wave::divideEPS (vector<INT16>canal) {

	vector<int>si;
	vector<double>coder;
	vector<double>decoder;
	vector<double>matrixX;
	vector<double>matrixP;
	vector<double>vectorEPS;
	vector<double>decoderEntropia;

	int b = 12;
	int k = 120 / r;
	int N = (ammount_of_samples / 2) -(k - 1) * ceil((ammount_of_samples / 2) / k);	//ceil(ammount_of_samples / 2)/k;
	double minLsr = 100;
	double Lsr;
	
	for (int p = 1; p <= k; p++) {
		double **A, *B, *X;
		int n, i, j;
		n = r;

		cout << setprecision(10) << fixed;
		A = new double *[n];
		B = new double[n];
		X = new double[n];

		for (i = 0; i < n; i++)
			A[i] = new double[n];

		double sumX = 0;
		double sumP = 0;

		for (i = 1; i <= r; i++) {
			for (j = 1; j <= r; j++) {
				for (int z = r + (N * p - N); z < N * p; z++) {
					sumX += canal.at(z - i) * canal.at(z - j);
					sumP += canal.at(z) * canal.at(z - i);
				}

				if (j == 1)
					matrixP.push_back(sumP);
				matrixX.push_back(sumX);
				sumX = 0;
				sumP = 0;
			}
		}

		int licznik = 0;
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				A[i][j] = matrixX.at(licznik);
				licznik++;
				B[i] = matrixP.at(i);
			}
		}

		if (ludist(n, A) && lusolve(n, A, B, X)) {}
		else cout << "DZIELNIK ZERO\n";

		for (int i = 0; i < r; i++) 
			vectorEPS.push_back(X[i]);

		for (i = 0; i < n; i++)
			delete[] A[i];
		delete[] A;
		delete[] B;
		delete[] X;

		auto max = max_element(begin(vectorEPS), end(vectorEPS));
		auto min = min_element(begin(vectorEPS), end(vectorEPS));
		auto position = 0;

		bool positive = true;
		if (abs(*min) > *max) {
			position = distance(begin(vectorEPS), min);
			*max = abs(*min);
			positive = false;
		}
		else {
			position = distance(begin(vectorEPS), max);
		}

		*max = float(*max);
		for (int i = 0; i < r; i++) {
			coder.push_back(floor(abs(vectorEPS.at(i)) / (*max) * (pow(2, b) - 1) + 0.5));
			si.push_back(sign(vectorEPS.at(i)));
		}

		for (int i = 0; i < r; i++) 
			decoder.push_back(((coder.at(i) / (pow(2, b) - 1)) * (*max)) * (si.at(i) * 2 - 1));
		
		decoderEntropia = sendEntropia(canal, decoder);
		Lsr = entro_minus(decoderEntropia) + ((32 + (r - 1) * (b + 1) + 10) / ammount_of_samples);

		if (minLsr > Lsr)
			minLsr = Lsr;

		si.clear();
		coder.clear();
		decoder.clear();
		matrixX.clear();
		matrixP.clear();
		vectorEPS.clear();
		decoderEntropia.clear();
	}

	return minLsr;
}

