#include "stdafx.h"
#include "WaveProcess.h"

int r = 7;
void start()
{
	WaveProcess wave1("ATrain.wav");
	/*WaveProcess wave2("BeautySlept.wav");
	WaveProcess wave3("death2.wav");
	WaveProcess wave4("experiencia.wav");
	WaveProcess wave5("chanchan.wav");
	WaveProcess wave6("female_speech.wav");
	WaveProcess wave7("FloorEssence.wav");
	WaveProcess wave8("ItCouldBeSweet.wav");
	WaveProcess wave9("Layla.wav");
	WaveProcess wave10("LifeShatters.wav");
	WaveProcess wave11("macabre.wav");
	WaveProcess wave12("male_speech.wav");
	WaveProcess wave13("SinceAlways.wav");
	WaveProcess wave14("thear1.wav");
	WaveProcess wave15("TomsDiner.wav");
	WaveProcess wave16("velvet.wav");

	vector <double> tab;
	double average = 0;

	tab.push_back(wave1.GetAverage());
	tab.push_back(wave2.GetAverage());
	tab.push_back(wave3.GetAverage());
	tab.push_back(wave4.GetAverage());
	tab.push_back(wave5.GetAverage());
	tab.push_back(wave6.GetAverage());
	tab.push_back(wave7.GetAverage());
	tab.push_back(wave8.GetAverage());
	tab.push_back(wave9.GetAverage());
	tab.push_back(wave10.GetAverage());
	tab.push_back(wave11.GetAverage());
	tab.push_back(wave12.GetAverage());
	tab.push_back(wave13.GetAverage());
	tab.push_back(wave14.GetAverage());
	tab.push_back(wave15.GetAverage());
	tab.push_back(wave16.GetAverage());

	for (int i = 0; i < 16; i++) {
		average += tab[i];
	}
	average /= 16;
	cout << average << endl;
	return average;*/
}
int main()
{
	start();
	//cout << "r=" << r << endl << endl;
	////for (r = 3; r < 10; r++)
	////{
	//	

	//	ofstream averageFile("average_test.txt", ios::app);
	//	averageFile << r << " " << start() << endl;
	////}

	system("pause");
	return 0;
}