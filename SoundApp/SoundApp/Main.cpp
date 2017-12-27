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

	double tab[32];
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

	tab[16] = wave1.getAverageEPS();
	tab[17] = wave2.getAverageEPS();
	tab[18] = wave3.getAverageEPS();
	tab[19] = wave4.getAverageEPS();
	tab[20] = wave5.getAverageEPS();
	tab[21] = wave6.getAverageEPS();
	tab[22] = wave7.getAverageEPS();
	tab[23] = wave8.getAverageEPS();
	tab[24] = wave9.getAverageEPS();
	tab[25] = wave10.getAverageEPS();
	tab[26] = wave11.getAverageEPS();
	tab[27] = wave12.getAverageEPS();
	tab[28] = wave13.getAverageEPS();
	tab[29] = wave14.getAverageEPS();
	tab[30] = wave15.getAverageEPS();
	tab[31] = wave16.getAverageEPS();

	for (int i = 0; i < 16; i++)
		average += tab[i];

	double average2 = 0;
	for (int i = 16; i < 32; i++)
		average2 += tab[i];

	average /= 16;
	average2 /= 16;
	ofstream averageFile("average.txt", ios::app);
	averageFile << r << " " << average << " " << average2 << endl;
}

void main() {

	/*for (int r = 2; r <= 40; r += 1) {
		cout << "r = " << r << endl;
		start(r);
	}*/

	Wave wave1("ATrain.wav", 3);

	/*start(3);
	start(4);
	start(5);
	start(6);
	start(8);
	start(10);
	start(12);
	start(15);
	start(20);
	start(24);
	start(30);
	start(40);
	start(60);
	start(120);*/
	
	system("pause");
}