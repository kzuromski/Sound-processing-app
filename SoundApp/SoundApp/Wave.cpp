#include "Wave.h"

Wave::Wave(string name_of_wave, int r) {
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
	ammount_of_samples = ((size_of_file - 36) - 2) / 2;

	NormalCanal();
	DifferentialCoder();

	//new_file << name_of_wave << endl;
	//new_file << "Liczba próbek z obu kana³ów: " << fixed << ammount_of_samples << endl;
	//double d = ((Energy(ammount_of_samples / 2, leftCanal) + Energy(ammount_of_samples / 2, rightCanal)) / 2);
	//new_file << "Przeciêtna energia sygna³u: " << fixed << d << endl; 
	//d = ((EnergyDifferential(ammount_of_samples / 2, leftDifferential) + EnergyDifferential(ammount_of_samples / 2, rightDifferential)) / 2);
	//new_file << "Przeciêtna energia sygna³u po skanowaniu ró¿nicowym: " << fixed << d << endl;  // drugie obliczenia
	//new_file << "Entropia: " << ((Entropy(leftCanal) + Entropy(rightCanal)) / 2) << endl; // Entropy dla normalnych
	//new_file << "Entropia z danych po skanowaniu ró¿nicowym: " << ((EntropyDifferential(leftDifferential) + EntropyDifferential(rightDifferential)) / 2) << endl;

	//r<2;40>
	//averageEPS = (SystemOfEquations(rightCanal) + SystemOfEquations(leftCanal)) / 2;
	//new_file << "Entropia ze wspó³czynnikiem " << averageEPS << endl;

	//r<10;600>
	//averageBit = (EntroBit(rightCanal) + EntroBit(leftCanal)) / 2;
	
	/*for (int i = 0; i < 2; i++)
		averageLsr += minLsrVector.at(i);
	averageLsr /= 2;*/

	//r<120;4>
	//averageLsr = (DivideEPS(rightCanal) + DivideEPS(leftCanal)) / 2;

	//DecoderDifferential(leftDifferential);
	//DecoderPredictive(leftCanal);

	BothCanals();
	averageLsr = (EntropyDifferential(predictCoderRight) + EntropyDifferential(predictCoderLeft)) / 2;
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

void Wave::ReadData() {
	fread(&riff, sizeof(FOURCC), 1, wf);
	fread(&size_of_file, sizeof(DWORD), 1, wf);
	fread(&wave, sizeof(FOURCC), 1, wf);
	fread(&fmt, sizeof(FOURCC), 1, wf);
	fread(&chunk, sizeof(FOURCC), 1, wf);
	fread(&pcm, sizeof(DWORD), 1, wf);
	fread(&chanel, sizeof(WORD), 1, wf);
	fread(&sample_rate, sizeof(DWORD), 1, wf);
	fread(&bytes_per_sec, sizeof(DWORD), 1, wf);
	fread(&block_alignment, sizeof(WORD), 1, wf);
	fread(&bits_per_sample, sizeof(WORD), 1, wf);
	fread(&data, sizeof(FOURCC), 1, wf);
	fread(&size_of_data, sizeof(DWORD), 1, wf);
}

void Wave::NormalCanal() {
	for (int i = 0; i < ammount_of_samples; i++){
		fread(&data_of_file, sizeof(INT16), 1, wf);
		if (i % 2 != 0)
			leftCanal.push_back(data_of_file);
		else
			rightCanal.push_back(data_of_file);
	}
}

void Wave::DifferentialCoder(){
	for (int i = 0; i < leftCanal.size(); i++){
		if (i == 0)
			leftDifferential.push_back(leftCanal.at(i));
		else
			leftDifferential.push_back(leftCanal.at(i) - leftCanal.at(i - 1));
	}

	for (int i = 0; i < rightCanal.size(); i++){
		if (i == 0)
			rightDifferential.push_back(rightCanal.at(i));
		else
			rightDifferential.push_back(rightCanal.at(i) - rightCanal.at(i - 1));
	}
}

vector<double> Wave::predictCoder(vector<INT16>canal, vector<double>vectorEPS) {
	double sumPredict = 0;
	vector <double> predictCoder;
	vector <double> predictValue;
	for (size_t i = 0; i < ammount_of_samples / 2; i++) {
		sumPredict = 0;
		if (i == 0)
			predictValue.push_back(canal.at(i));
		else if (i < r)
			predictValue.push_back(canal.at(i) - canal.at(i - 1));
		else {
			for (size_t j = 1; j <= r; j++)
				sumPredict += vectorEPS.at(j - 1) * canal.at(i - j);
			if (sumPredict > 32768 - 1)
				sumPredict = 32768 - 1;
			else if (sumPredict < -32768)
				sumPredict = -32768;
			predictValue.push_back(floor(sumPredict + 0.5));
		}
	}

	for (int i = 0; i < ammount_of_samples / 2; i++) {
		if (i < r)
			predictCoder.push_back(predictValue.at(i));
		else
			predictCoder.push_back(canal.at(i) - predictValue.at(i));

	}

	return predictCoder;
}

double Wave::Energy(double a, vector<INT16> b) {
	double full = 0;
	for (INT32 i = 0; i < a; i++)
		full = (double)(full + (((double)b.at(i) * (double)b.at(i))));

	full = full / a;

	return full;
}

double Wave::EnergyDifferential(double a, vector<double> b) {
	double full = 0;
	for (INT32 i = 0; i < a; i++)
		full = (double)(full + (b.at(i) * b.at(i)));

	full = full / a;

	return full;
}

double Wave::Entropy(vector<INT16> a) {
	double Entropy = 0;
	vector <double> buffor(131072, 0); 
	for (int i = 0; i < a.size(); i++)
		buffor.at(a.at(i) + 65536)++; 

	for (int i = 0; i < 131072; i++){
		if (buffor.at(i) != 0) {
			double p_i = (double)buffor.at(i) / a.size();
			Entropy = Entropy + (p_i*log2(p_i));
		}
	}
	return Entropy *(-1);
}

double Wave::EntropyDifferential(vector<double> a){
	double Entropy = 0;
	vector <double> buffor(262144, 0); //131072
	for (int i = 0; i < a.size(); i++)
		buffor.at(a.at(i) + 131072)++; //65536
	for (int i = 0; i < 262144; i++){ //131072
		if (buffor.at(i) != 0){
			double p_i = (double)buffor.at(i) / a.size();
			Entropy = Entropy + (p_i*log2(p_i));
		}
	}
	return Entropy *(-1);
}

bool Wave::ludist(int n, double ** A){
	int i, j, k;

	for (k = 0; k < n - 1; k++){
		if (fabs(A[k][k]) < eps) return false;

		for (i = k + 1; i < n; i++)
			A[i][k] /= A[k][k];

		for (i = k + 1; i < n; i++)
			for (j = k + 1; j < n; j++)
				A[i][j] -= A[i][k] * A[k][j];
	}

	return true;
}

bool Wave::lusolve(int n, double ** A, double * B, double * X){
	int    i, j;
	double s;

	X[0] = B[0];

	for (i = 1; i < n; i++){
		s = 0;
		for (j = 0; j < i; j++) s += A[i][j] * X[j];
		X[i] = B[i] - s;
	}

	if (fabs(A[n - 1][n - 1]) < eps) return false;
	X[n - 1] /= A[n - 1][n - 1];
	for (i = n - 2; i >= 0; i--){
		s = 0;
		for (j = i + 1; j < n; j++) s += A[i][j] * X[j];
		if (fabs(A[i][i]) < eps) return false;
		X[i] = (X[i] - s) / A[i][i];
	}

	return true;
}

double Wave::SystemOfEquations(vector<INT16>canal) {
	double **A, *B, *X;
	int n = r;

	cout << setprecision(10) << fixed;
	A = new double *[n];
	B = new double[n];
	X = new double[n];

	for (int i = 0; i < n; i++)
		A[i] = new double[n];

	int N = ammount_of_samples / 2;
	double sumX = 0;
	double sumP = 0;
	vector<double>matrixX;
	vector<double>matrixP;
	vector<double> vectorEPS;

	for (int i = 1; i <= r; i++) {
		for (int j = 1; j <= r; j++) {
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

	int counterVector = 0;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			A[i][j] = matrixX.at(counterVector);
			counterVector++;
			B[i] = matrixP.at(i);
		}
	}

	if (ludist(n, A) && lusolve(n, A, B, X)) {}
	else cout << "DZIELNIK ZERO\n";

	for (int i = 0; i < r; i++) 
		vectorEPS.push_back(X[i]);

	vector<double> codingSamples;
	codingSamples = predictCoder(canal, vectorEPS);
	double returnEntropia = EntropyDifferential(codingSamples);

	for (int i = 0; i < n; i++)
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

double Wave::EntroBit(vector<INT16>canal) {
	double **A, *B, *X;
	int n = r;
	int N = ammount_of_samples / 2;

	cout << setprecision(10) << fixed;
	A = new double *[n];
	B = new double[n];
	X = new double[n];

	for (int i = 0; i < n; i++)
		A[i] = new double[n];

	double sumX = 0;
	double sumP = 0;
	vector<double>matrixX;
	vector<double>matrixP;
	vector<double>vectorEPS;

	for (int i = 1; i <= r; i++) {
		for (int j = 1; j <= r; j++) {
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

	int counterVector = 0;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			A[i][j] = matrixX.at(counterVector);
			counterVector++;
			B[i] = matrixP.at(i);
		}
	}

	if (ludist(n, A) && lusolve(n, A, B, X)) {}
	else cout << "DZIELNIK ZERO\n";

	for (int i = 0; i < r; i++) 
		vectorEPS.push_back(X[i]);

	for (int i = 0; i < n; i++)
		delete[] A[i];
	delete[] A;
	delete[] B;
	delete[] X;

	double max = *max_element(begin(vectorEPS), end(vectorEPS));
	double min = *min_element(begin(vectorEPS), end(vectorEPS));

	if (abs(min) > max) 
		max = abs(min);
	max = float(max);

	vector<int>si;
	vector<double>scale;
	vector<double>descale;
	vector<double> codingSamples;

	double Lsr;
	double entropia = 0;
	double minLsr = 100;
	int diagramBit = 0;
	for (int b = 5; b <= 16; b++) {

		for (int i = 0; i < r; i++) {
			scale.push_back(floor(abs(vectorEPS.at(i)) / (max) * (pow(2, b) - 1) + 0.5));
			si.push_back(sign(vectorEPS.at(i)));
		}

		for (int i = 0; i < r; i++)
			descale.push_back(((scale.at(i) / (pow(2, b) - 1)) * (max)) * (si.at(i) * 2 - 1));

		codingSamples = predictCoder(canal, descale);
		Lsr = EntropyDifferential(codingSamples) + ((32 + (r - 1) * (b + 1) + 10) / ammount_of_samples);

		if (minLsr > Lsr) {
			minLsr = Lsr;
			diagramBit = b;
		}	

		si.clear();
		scale.clear();
		descale.clear();
		codingSamples.clear();
	}
	minLsrVector.push_back(minLsr);
	
	return diagramBit;
}

double Wave::DivideEPS (vector<INT16>canal) {
	vector<int>si;
	vector<double>scale;
	vector<double>descale;
	vector<double>matrixX;
	vector<double>matrixP;
	vector<double>vectorEPS;
	vector<double>codingSamples;

	int b = 12;
	int k = 120 / r;
	double N = (ammount_of_samples / 2) / k;//(ammount_of_samples / 2) - (k - 1) * ceil((ammount_of_samples / 2) / k);
	double minLsr = 100;
	double Lsr;
	
	for (int p = 1; p <= k; p++) {
		double **A, *B, *X;
		int n = r;

		cout << setprecision(10) << fixed;
		A = new double *[n];
		B = new double[n];
		X = new double[n];

		for (int i = 0; i < n; i++)
			A[i] = new double[n];

		double sumX = 0;
		double sumP = 0;

		for (int i = 1; i <= r; i++) {
			for (int j = 1; j <= r; j++) {
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

		int counterVector = 0;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				A[i][j] = matrixX.at(counterVector);
				counterVector++;
				B[i] = matrixP.at(i);
			}
		}

		if (ludist(n, A) && lusolve(n, A, B, X)) {}
		else cout << "DZIELNIK ZERO\n";

		for (int i = 0; i < r; i++)
			vectorEPS.push_back(X[i]);

		for (int i = 0; i < n; i++)
			delete[] A[i];
		delete[] A;
		delete[] B;
		delete[] X;

		double max = *max_element(begin(vectorEPS), end(vectorEPS));
		double min = *min_element(begin(vectorEPS), end(vectorEPS));

		if (abs(min) > max) 
			max = abs(min);
		max = float(max);

		for (int i = 0; i < r; i++) {
			scale.push_back(floor(abs(vectorEPS.at(i)) / (max) * (pow(2, b) - 1) + 0.5));
			si.push_back(sign(vectorEPS.at(i)));
		}

		for (int i = 0; i < r; i++) 
			descale.push_back(((scale.at(i) / (pow(2, b) - 1)) * (max)) * (si.at(i) * 2 - 1));

		codingSamples = predictCoder(canal, descale);
		Lsr = EntropyDifferential(codingSamples) + ((32 + (r - 1) * (b + 1) + 10) / ammount_of_samples);

		if (minLsr > Lsr)
			minLsr = Lsr;

		si.clear();
		scale.clear();
		descale.clear();
		vectorEPS.clear();
		codingSamples.clear();
	}

	return minLsr;
}

void Wave::DecoderDifferential(vector<double>canal) {
	vector <double> decoderCanal;
	for (int i = 0; i < canal.size(); i++){
		if (i == 0)
			decoderCanal.push_back(canal.at(i));
		else
			decoderCanal.push_back(canal.at(i) + decoderCanal.at(i - 1));
	}
}

void Wave::DecoderPredictive(vector<INT16>canal) {
	double **A, *B, *X;
	int n = r;

	cout << setprecision(5) << fixed;
	A = new double *[n];
	B = new double[n];
	X = new double[n];

	for (int i = 0; i < n; i++)
		A[i] = new double[n];

	int N = ammount_of_samples / 2;
	double sumX = 0;
	double sumP = 0;
	vector<double>matrixX;
	vector<double>matrixP;
	vector<double> vectorEPS;

	for (int i = 1; i <= r; i++) {
		for (int j = 1; j <= r; j++) {
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

	int counterVector = 0;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			A[i][j] = matrixX.at(counterVector);
			counterVector++;
			B[i] = matrixP.at(i);
		}
	}

	if (ludist(n, A) && lusolve(n, A, B, X)) {}
	else cout << "DZIELNIK ZERO\n";

	for (int i = 0; i < r; i++)
		vectorEPS.push_back(X[i]);

	for (int i = 0; i < n; i++)
		delete[] A[i];
	delete[] A;
	delete[] B;
	delete[] X;

	vector <double> coderCanal = predictCoder(canal, vectorEPS);
	vector <double> decoderCanal;
	double sumPredict = 0;
	vector <double> predictValue;

	for (size_t i = 0; i < ammount_of_samples / 2; i++) {
		sumPredict = 0;
		if (i == 0)
			predictValue.push_back(canal.at(i));
		else if (i < r)
			predictValue.push_back(canal.at(i) - canal.at(i - 1));
		else {
			for (size_t j = 1; j <= r; j++)
				sumPredict += vectorEPS.at(j - 1) * canal.at(i - j);
			if (sumPredict > 32768 - 1)
				sumPredict = 32768 - 1;
			else if (sumPredict < -32768)
				sumPredict = -32768;
			predictValue.push_back(floor(sumPredict + 0.5));
		}
	}

	for (int i = 0; i < ammount_of_samples / 2; i++) {
		if (i == 0)
			decoderCanal.push_back(predictValue.at(i));
		else if (i < r)
			decoderCanal.push_back(predictValue.at(i) + canal.at(i - 1));
		else
			decoderCanal.push_back(coderCanal.at(i) + predictValue.at(i));
	}
}

void Wave::BothCanals() {

	int r1 = r * 0.5; // 0.5 <=> 0.67
	int r2 = r * 0.5; // 0.5 <=> 0.34

	//RIGHT---------------------------------------------------------------------------

	double **A1, *B1, *X1;
	int n1 = r;

	cout << setprecision(5) << fixed;
	A1 = new double *[n1];
	B1 = new double[n1];
	X1 = new double[n1];

	for (int i = 0; i < n1; i++)
		A1[i] = new double[n1];

	int N1 = ammount_of_samples / 2;
	double sumX1 = 0;
	double sumP1 = 0;
	vector<double>matrixXfirst;
	vector<double>matrixPfirst;
	vector<double> vectorEPSfirst;

	for (int i = 1; i <= r; i++) {
		for (int j = 1; j <= r; j++) {
			for (int z = r; z < N1; z++) {
				if(i > r1 && j > r1)
					sumX1 += leftCanal.at(z - i + r1) * leftCanal.at(z - j + r1);
				else if (i > r1)
					sumX1 += leftCanal.at(z - i + r1) * rightCanal.at(z - j);
				else if (j > r1) 
					sumX1 += rightCanal.at(z - i) * leftCanal.at(z - j + r1);
				else
					sumX1 += rightCanal.at(z - i) * rightCanal.at(z - j);
				sumP1 += rightCanal.at(z) * rightCanal.at(z - i);
			}

			if (j == 1)
				matrixPfirst.push_back(sumP1);
			matrixXfirst.push_back(sumX1);
			sumX1 = 0;
			sumP1 = 0;
		}
	}

	int counterVector1 = 0;
	for (int i = 0; i < n1; i++) {
		for (int j = 0; j < n1; j++) {
			A1[i][j] = matrixXfirst.at(counterVector1);
			counterVector1++;
			B1[i] = matrixPfirst.at(i);
		}
	}

	if (ludist(n1, A1) && lusolve(n1, A1, B1, X1)) {}
	else cout << "DZIELNIK ZERO\n";

	for (int i = 0; i < r; i++)
		vectorEPSfirst.push_back(X1[i]);

	for (int i = 0; i < n1; i++)
		delete[] A1[i];
	delete[] A1;
	delete[] B1;
	delete[] X1;

	//LEFT---------------------------------------------------------------------------

	double **A, *B, *X;
	int n = r;

	cout << setprecision(5) << fixed;
	A = new double *[n];
	B = new double[n];
	X = new double[n];

	for (int i = 0; i < n; i++)
		A[i] = new double[n];

	int N = ammount_of_samples / 2;
	double sumX = 0;
	double sumP = 0;
	vector<double>matrixXsecond;
	vector<double>matrixPseconds;
	vector<double>vectorEPSsecond;

	for (int i = 1; i <= r; i++) {
		for (int j = 1; j <= r; j++) {
			for (int z = r; z < N; z++) {
				if (i > r2 && j > r2)
					sumX1 += leftCanal.at(z - i + r2 - 1) * leftCanal.at(z - j + r2 - 1);
				else if (i > r2) 
					sumX += leftCanal.at(z - i + r2 - 1) * rightCanal.at(z - j);
				else if (j > r2) 
					sumX += rightCanal.at(z - i) * leftCanal.at(z - j + r2 - 1);
				else
					sumX += rightCanal.at(z - i) * rightCanal.at(z - j);
				sumP += rightCanal.at(z) * rightCanal.at(z - i);
			}

			if (j == 1)
				matrixPseconds.push_back(sumP);
			matrixXsecond.push_back(sumX);
			sumX = 0;
			sumP = 0;
		}
	}

	int counterVector = 0;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			A[i][j] = matrixXsecond.at(counterVector);
			counterVector++;
			B[i] = matrixPseconds.at(i);
		}
	}

	if (ludist(n, A) && lusolve(n, A, B, X)) {}
	else cout << "DZIELNIK ZERO\n";

	for (int i = 0; i < r; i++)
		vectorEPSsecond.push_back(X[i]);

	for (int i = 0; i < n; i++)
		delete[] A[i];
	delete[] A;
	delete[] B;
	delete[] X;
	
	//RIGHT-------------------------------------------------------------------------- -

	double sumPredict1 = 0;
	vector <double> predictCoderRightfirst;
	vector <double> predictValueRightfirst;
	for (size_t i = 0; i < ammount_of_samples / 2; i++) {
		sumPredict1 = 0;
		if (i == 0)
			predictValueRightfirst.push_back(rightCanal.at(i));
		else if (i < r)
			predictValueRightfirst.push_back(rightCanal.at(i) - rightCanal.at(i - 1));
		else {
			for (size_t j = 1; j <= r1; j++)
				sumPredict1 += vectorEPSfirst.at(j - 1) * rightCanal.at(i - j);
			if (sumPredict1 > 32768 - 1)
				sumPredict1 = 32768 - 1;
			else if (sumPredict1 < -32768)
				sumPredict1 = -32768;
			predictValueRightfirst.push_back(floor(sumPredict1 + 0.5));
		}
	}

	for (int i = 0; i < ammount_of_samples / 2; i++) {
		if (i < r)
			predictCoderRightfirst.push_back(predictValueRightfirst.at(i));
		else
			predictCoderRightfirst.push_back(rightCanal.at(i) - predictValueRightfirst.at(i));

	}

	//RIGHT-------------------------------------------------------------------------- -

	double sumPredict2 = 0;
	vector <double> predictCoderRightsecond;
	vector <double> predictValueRightsecond;
	for (size_t i = 0; i < ammount_of_samples / 2; i++) {
		sumPredict2 = 0;
		if (i == 0)
			predictValueRightsecond.push_back(rightCanal.at(i));
		else if (i < r)
			predictValueRightsecond.push_back(rightCanal.at(i) - rightCanal.at(i - 1));
		else {
			for (size_t j = 0; j <= r2 - 1; j++)
			sumPredict2 += vectorEPSsecond.at(j) * rightCanal.at(i - j);
			if (sumPredict2 > 32768 - 1)
				sumPredict2 = 32768 - 1;
			else if (sumPredict2 < -32768)
				sumPredict2 = -32768;
			predictValueRightsecond.push_back(floor(sumPredict2 + 0.5));
		}
	}

	for (int i = 0; i < ammount_of_samples / 2; i++) {
		if (i < r)
			predictCoderRightsecond.push_back(predictValueRightsecond.at(i));
		else
			predictCoderRightsecond.push_back(rightCanal.at(i) - predictValueRightsecond.at(i));

	}

	//LEFT---------------------------------------------------------------------------

	double sumPredict = 0;
	vector <double> predictCoderLeftfirst;
	vector <double> predictValueLeftfirst;
	for (size_t i = 0; i < ammount_of_samples / 2; i++) {
		sumPredict = 0;
		if (i == 0)
			predictValueLeftfirst.push_back(leftCanal.at(i));
		else if (i < r)
			predictValueLeftfirst.push_back(leftCanal.at(i) - leftCanal.at(i - 1));
		else {
			for (size_t j = 1; j <= r2; j++)
				sumPredict += vectorEPSsecond.at(j - 1) * leftCanal.at(i - j);
			if (sumPredict > 32768 - 1)
				sumPredict = 32768 - 1;
			else if (sumPredict < -32768)
				sumPredict = -32768;
			predictValueLeftfirst.push_back(floor(sumPredict + 0.5));
		}
	}

	for (int i = 0; i < ammount_of_samples / 2; i++) {
		if (i < r)
			predictCoderLeftfirst.push_back(predictValueLeftfirst.at(i));
		else
			predictCoderLeftfirst.push_back(leftCanal.at(i) - predictValueLeftfirst.at(i));

	}

	//LEFT---------------------------------------------------------------------------

	double sumPredict3 = 0;
	vector <double> predictCoderLeftsecond;
	vector <double> predictValueLeftsecond;
	for (size_t i = 0; i < ammount_of_samples / 2; i++) {
		sumPredict3 = 0;
		if (i == 0)
			predictValueLeftsecond.push_back(leftCanal.at(i));
		else if (i < r)
			predictValueLeftsecond.push_back(leftCanal.at(i) - leftCanal.at(i - 1));
		else {
			for (size_t j = 1; j <= r1; j++)
				sumPredict3 += vectorEPSfirst.at(j - 1) * leftCanal.at(i - j);
			if (sumPredict3 > 32768 - 1)
				sumPredict3 = 32768 - 1;
			else if (sumPredict3 < -32768)
				sumPredict3 = -32768;
			predictValueLeftsecond.push_back(floor(sumPredict3 + 0.5));
		}
	}

	for (int i = 0; i < ammount_of_samples / 2; i++) {
		if (i < r)
			predictCoderLeftsecond.push_back(predictValueLeftsecond.at(i));
		else
			predictCoderLeftsecond.push_back(leftCanal.at(i) - predictValueLeftsecond.at(i));

	}
	
	//double sum = 0;
	//for (int i = 0; i < 6; i++) {
	//	cout << vectorEPSfirst.at(i) << endl;
	//	cout << vectorEPSsecond.at(i) << endl;
	//	sum += vectorEPSsecond.at(i);
	//}
	//cout << sum << endl << endl;

	for (int i = 0; i < ammount_of_samples / 2; i++) {
		predictCoderLeft.push_back(predictCoderRightfirst.at(i) + predictCoderLeftfirst.at(i));
		predictCoderRight.push_back(predictCoderLeftsecond.at(i) + predictCoderRightsecond.at(i));
	}
}