#include "Wave.h"

void start(int r) {
	Wave wave1("ATrain.wav", r);
	Wave wave2("BeautySlept.wav", r);
	Wave wave3("death2.wav", r);
	Wave wave4("experiencia.wav", r);
	Wave wave5("chanchan.wav", r);
	Wave wave6("female_speech.wav", r);
	Wave wave7("FloorEssence.wav", r);
	Wave wave8("ItCouldBeSweet.wav", r);
	Wave wave9("Layla.wav", r);
	Wave wave10("LifeShatters.wav", r);
	Wave wave11("macabre.wav", r);
	Wave wave12("male_speech.wav", r);
	Wave wave13("SinceAlways.wav", r);
	Wave wave14("thear1.wav", r);
	Wave wave15("TomsDiner.wav", r);
	Wave wave16("velvet.wav", r);

	double tab[16];
	double average = 0;

	tab[0] = wave1.getAverageBit();
	tab[1] = wave2.getAverageBit();
	tab[2] = wave3.getAverageBit();
	tab[3] = wave4.getAverageBit();
	tab[4] = wave5.getAverageBit();
	tab[5] = wave6.getAverageBit();
	tab[6] = wave7.getAverageBit();
	tab[7] = wave8.getAverageBit();
	tab[8] = wave9.getAverageBit();
	tab[9] = wave10.getAverageBit();
	tab[10] = wave11.getAverageBit();
	tab[11] = wave12.getAverageBit();
	tab[12] = wave13.getAverageBit();
	tab[13] = wave14.getAverageBit();
	tab[14] = wave15.getAverageBit();
	tab[15] = wave16.getAverageBit();

	for (int i = 0; i < 16; i++)
		average += tab[i];

	average /= 16;
	ofstream averageFile("average.txt", ios::app);
	averageFile << r << " " << average << endl;
}

void main() {

	for (int r = 10; r <= 600; r += 10) {
		cout << "r = " << r << endl;
		start(r);
	}

	system("pause");
}